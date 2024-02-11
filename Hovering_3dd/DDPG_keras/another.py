import gym
import tensorflow as tf
from tensorflow.keras import layers
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pickle
import csv
from mpl_toolkits import mplot3d

import glob
import ctypes
import math
### Setup the environment

num_states=46 # 9errors, 3dt, 4p'_r'_p''_r'', 30pr_inputs

num_actions=3 # u_x, u_y, u_z

mass=0.67 #in kg

prev_input_number = 30

g=9.81

upper_bound= mass*g # i.e. coefficient of thrust to weight=2

lower_bound= -mass * g # we cannot do better without upside-down

vel_reward_weight_factor=0

libfile = glob.glob('../Sim_intfc/out/build/linux-default/libinterface_simulink.so.1')[0]
mylib = ctypes.CDLL(libfile)

class EnvUAV_cmd_z:
    simulink_input_length=3

    simulink_states_length=37
    print_period=1 #in seconds
    episode_timeout=3 #in seconds
    termination_causes={'passed window plane':0,'time_out':1,'nan_error':2}
    def __init__(self,simulation_dt,sampling_dt,mass=1,CTW=2,g=9.81):
        mylib.LoadModel()
        mylib.run_a_step.restype = np.ctypeslib.ndpointer(dtype=np.float64,shape=(self.simulink_states_length,))
        mylib.run_a_step.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64),ctypes.c_int]

        self.simulation_dt=simulation_dt
        self.sampling_dt=sampling_dt
        self.sim_steps_per_samp_period=round(self.sampling_dt/self.simulation_dt) # Note sim dt must not be larger than sampling dt

        self.steps_per_print_period=int(EnvUAV_cmd_z.print_period/self.sampling_dt)
        self.sampling_step_counter=0
        self.reset()
    
    def __call__(self,input):
        self.sampling_step_counter=self.sampling_step_counter+1
        inputs = np.zeros(EnvUAV_cmd_z.simulink_input_length,dtype=np.float64)
        inputs[:3] = np.copy(input[:3])
        
        for i in range(self.sim_steps_per_samp_period): # Simulink solver is run at a higher rate for numerical stability
            self.states_meas = mylib.run_a_step(inputs,len(inputs))
        
        for i in range(9): #errors
            self.all_states_meas[i] = np.copy(self.states_meas[i] - self.states_meas[i + 28])
        
        indices_to_copy = [25, 26, 27, 16, 17, 19, 20]
        for i, index in enumerate(indices_to_copy, start=9):# dt,omegas,omegas'
            self.all_states_meas[i] = np.copy(self.states_meas[index])

        if num_states==46: # d =10 for each input (30)
            self.all_states_meas[16:19] = np.copy(input[:3])
            self.all_states_meas[19:46] = np.copy(self.prevv_input)
            self.prevv_input[3:27] = self.prevv_input[:24]
            self.prevv_input[:3] = np.copy(input[:3])

        #cost function
        elements = np.array([
            self.all_states_meas[0],
            self.all_states_meas[1],
            self.all_states_meas[2],
            vel_reward_weight_factor * self.all_states_meas[3],
            vel_reward_weight_factor * self.all_states_meas[4],
            vel_reward_weight_factor * self.all_states_meas[5]
        ])
        cost = np.linalg.norm(elements)       
        self.reward=-cost*self.sampling_dt*100 # By 100 to 25 up reward for training
        
        if en_goal_termination==True:    
            if (self.all_states_meas[0]< 0):
                self.termination_cause=EnvUAV_cmd_z.termination_causes['passed window plane']
                print('Passed window plane!')
                self.done=True
            
        if (self.sampling_step_counter*self.sampling_dt>EnvUAV_cmd_z.episode_timeout):
            self.termination_cause=EnvUAV_cmd_z.termination_causes['time_out']
            self.done=True
        
        if any(math.isnan(self.all_states_meas[i]) for i in range(6)):            
            self.all_states_meas=np.zeros(num_states)
            self.reward=np.ones(1)*0.0
            self.termination_cause=EnvUAV_cmd_z.termination_causes['nan_error']
            print('Warning NaN termination')
            self.done=True
        return self.all_states_meas, self.reward, self.done, self.termination_cause, self.states_meas[28:31], self.states_meas[0:3]

    def reset(self):
        if self.sampling_step_counter>0:
            mylib.TermModel()
        self.sampling_step_counter=0
        mylib.LoadModel()
        mylib.InitModel()
        self.done=False
        self.all_states_meas=np.zeros(num_states)
        self.prevv_input=np.zeros(prev_input_number-3)
        self.termination_cause=EnvUAV_cmd_z.termination_causes['passed window plane']
        return self.all_states_meas

    def print_states_periodic(self):
        if self.sampling_step_counter % self.steps_per_print_period ==1: #here
            print('simulation time (s): ' + str(self.sampling_step_counter*self.sampling_dt))
            print('===================')

simulink_dt=0.001 # Simulink solver step time
sampling_dt=0.005 # The one that would be used experimentally. e.g. 0.005 corresponds to 200Hz

UAV_alt_ctrl=EnvUAV_cmd_z(simulink_dt,sampling_dt)


print(UAV_alt_ctrl.sim_steps_per_samp_period)

class OUActionNoise:
    def __init__(self, mean, std_deviation, theta=0.15, dt=1e-2, x_initial=None):
        self.theta = theta
        self.mean = mean
        self.std_dev = std_deviation
        self.dt = dt
        self.x_initial = x_initial
        self.reset()

    def __call__(self):
        # Formula taken from https://www.wikipedia.org/wiki/Ornstein-Uhlenbeck_process.
        x = (
            self.x_prev
            + self.theta * (self.mean - self.x_prev) * self.dt
            + self.std_dev * np.sqrt(self.dt) * np.random.normal(size=self.mean.shape)
        )
        # Store x into x_prev
        # Makes next noise dependent on current one
        self.x_prev = x
        return x

    def reset(self):
        if self.x_initial is not None:
            self.x_prev = self.x_initial
        else:
            self.x_prev = np.zeros_like(self.mean)
            
class Buffer:
    def __init__(self, buffer_capacity=100000, batch_size=64):
        # Number of "experiences" to store at max
        self.buffer_capacity = buffer_capacity
        # Num of tuples to train on.
        self.batch_size = batch_size

        # Its tells us num of times record() was called.
        self.buffer_counter = 0

        # Instead of list of tuples as the exp.replay concept go
        # We use different np.arrays for each tuple element
        self.state_buffer = np.zeros((self.buffer_capacity, num_states))
        self.action_buffer = np.zeros((self.buffer_capacity, num_actions))
        self.reward_buffer = np.zeros((self.buffer_capacity, 1)) 
        self.next_state_buffer = np.zeros((self.buffer_capacity, num_states))

    # Takes (s,a,r,s') obervation tuple as input
    def record(self, obs_tuple):
        # Set index to zero if buffer_capacity is exceeded,
        # replacing old records
        index = self.buffer_counter % self.buffer_capacity

        self.state_buffer[index] = obs_tuple[0]
        self.action_buffer[index] = obs_tuple[1]
        # self.action_buffer[index] = np.array(obs_tuple[1]).flatten()
        self.reward_buffer[index] = obs_tuple[2]
        self.next_state_buffer[index] = obs_tuple[3]

        self.buffer_counter += 1

    # Eager execution is turned on by default in TensorFlow 2. Decorating with tf.function allows
    # TensorFlow to build a static graph out of the logic and computations in our function.
    # This provides a large speed up for blocks of code that contain many small TensorFlow operations such as this one.
    @tf.function
    def update(
        self, state_batch, action_batch, reward_batch, next_state_batch,
    ):
        # Training and updating Actor & Critic networks.
        # See Pseudo Code.
        with tf.GradientTape() as tape:
            target_actions = target_actor(next_state_batch, training=True)
            y = reward_batch + gamma * target_critic(
                [next_state_batch, target_actions], training=True
            )
            critic_value = critic_model([state_batch, action_batch], training=True)
            critic_loss = tf.math.reduce_mean(tf.math.square(y - critic_value))

        critic_grad = tape.gradient(critic_loss, critic_model.trainable_variables)
        critic_optimizer.apply_gradients(
            zip(critic_grad, critic_model.trainable_variables)
        )

        with tf.GradientTape() as tape:
            actions = actor_model(state_batch, training=True)
            critic_value = critic_model([state_batch, actions], training=True)
            # Used `-value` as we want to maximize the value given
            # by the critic for our actions
            actor_loss = -tf.math.reduce_mean(critic_value)

        actor_grad = tape.gradient(actor_loss, actor_model.trainable_variables)
        actor_optimizer.apply_gradients(
            zip(actor_grad, actor_model.trainable_variables)
        )

    # We compute the loss and update parameters
    def learn(self):
        # Get sampling range
        record_range = min(self.buffer_counter, self.buffer_capacity)
        # Randomly sample indices
        batch_indices = np.random.choice(record_range, self.batch_size)

        # Convert to tensors
        state_batch = tf.convert_to_tensor(self.state_buffer[batch_indices])
        action_batch = tf.convert_to_tensor(self.action_buffer[batch_indices])
        reward_batch = tf.convert_to_tensor(self.reward_buffer[batch_indices])
        reward_batch = tf.cast(reward_batch, dtype=tf.float32)
        next_state_batch = tf.convert_to_tensor(self.next_state_buffer[batch_indices])

        self.update(state_batch, action_batch, reward_batch, next_state_batch)


# This update target parameters slowly
# Based on rate `tau`, which is much less than one.
@tf.function
def update_target(target_weights, weights, tau):
    for (a, b) in zip(target_weights, weights):
        a.assign(b * tau + a * (1 - tau))

def get_actor():
    # Initialize weights between -3e-3 and 3-e3
    last_init = tf.random_uniform_initializer(minval=-0.2, maxval=0.2) #MC
    multiple_factor=2
    inputs = layers.Input(shape=(num_states,))
    out = layers.Dense(multiple_factor*256, activation="relu")(inputs)
    out = layers.Dense(multiple_factor*256, activation="relu")(out)
    outputs = layers.Dense(3, activation="tanh", kernel_initializer=last_init)(out)

    # Our upper bound is 2.0 for Pendulum.
    outputs = outputs * upper_bound
    model = tf.keras.Model(inputs, outputs)
    return model


def get_critic():
    # State as input
    multiple_factor=2
    state_input = layers.Input(shape=(num_states))
    state_out = layers.Dense(multiple_factor*16, activation="relu")(state_input)
    state_out = layers.Dense(multiple_factor*32, activation="relu")(state_out)

    # Action as input
    action_input = layers.Input(shape=(num_actions))
    action_out = layers.Dense(multiple_factor*32, activation="relu")(action_input)

    # Both are passed through seperate layer before concatenating
    concat = layers.Concatenate()([state_out, action_out])

    out = layers.Dense(multiple_factor*256, activation="relu")(concat)
    out = layers.Dense(multiple_factor*256, activation="relu")(out)
    outputs = layers.Dense(1)(out)

    # Outputs single value for give state-action
    model = tf.keras.Model([state_input, action_input], outputs)

    return model

def policy(state, noise_object):
    sampled_actions = tf.squeeze(actor_model(state))
    noise = noise_object()
    # Adding noise to action
    sampled_actions = sampled_actions.numpy() + noise # WARNING NOISE MIGHT BE DISABLED

    # We make sure action is within bounds
    legal_action = np.clip(sampled_actions, lower_bound, upper_bound)

    return [np.squeeze(legal_action)]

# RL models
actor_model = get_actor()
critic_model = get_critic()

target_actor = get_actor()
target_critic = get_critic()

# Making the weights equal initially
target_actor.set_weights(actor_model.get_weights())
target_critic.set_weights(critic_model.get_weights())

# Learning rate for actor-critic models
critic_lr = 0.002
actor_lr = 0.001

critic_optimizer = tf.keras.optimizers.Adam(critic_lr)
actor_optimizer = tf.keras.optimizers.Adam(actor_lr)

total_episodes = 2
# Discount factor for future rewards
gamma = 0.99
# Used to update target networks
tau = 0.005

buffer = Buffer(400000, 1024)

# nosie model
init_std_dev = 0.4
final_exponent_val=10 #formula std_dev=init_std_dev * exp^(-(ep/total_episodes)*final_exponent_val)
ou_noise = OUActionNoise(mean=np.zeros(1), std_deviation=float(init_std_dev) * np.ones(1))
en_noise_decay=False

# model behaviour
en_goal_termination=False

# To store reward history of each episode
ep_reward_list = []

# To store average reward history of last few episodes
avg_reward_list = []

# To store final states of the system
final_states_reward_list = []  # [states, reward]

# To store termiantion status
termination_causes_list=[]

# To store all state_action pairs
states_action_list=[]
states_action_lists=[]
window_state_list=[]
window_state_lists=[]
drone_state_list=[]
drone_state_lists=[]

learning_frequency=256 #MC

for ep in range(total_episodes):

    #prev_state = env.reset()
    prev_state = UAV_alt_ctrl.reset()
    episodic_reward = 0 # Total rewards in one episode

    while True:
        # Uncomment this to see the Actor in action
        # But not in a python notebook.
        # env.render()

        

        tf_prev_state = tf.expand_dims(tf.convert_to_tensor(prev_state), 0) # Necessary for tensorflow (check immutable/mutable data types if you want to know more)

        action = policy(tf_prev_state, ou_noise)
        action = list(action[0]) 
        # Recieve state and reward from environment.
        #state, reward, done, info = env.step(action)
        state, reward, done, termination_cause, window_position, drone_position=UAV_alt_ctrl(action)
        states_action_list.append([state.copy(),action.copy()])
        window_state_list.append(window_position.copy())
        drone_state_list.append(drone_position.copy())

        buffer.record((prev_state.copy(), action.copy(), reward.copy(), state.copy()))
        episodic_reward += reward

        if buffer.buffer_counter>buffer.batch_size:
            if (buffer.buffer_counter%(int(buffer.batch_size/learning_frequency))==1):
                buffer.learn()
                update_target(target_actor.variables, actor_model.variables, tau)
                update_target(target_critic.variables, critic_model.variables, tau)


        # End this episode when `done` is True
        if done:
            break

        prev_state = np.copy(state)

    ep_reward_list.append(episodic_reward)
    termination_causes_list.append(termination_cause)
    final_states_reward_list.append([state.copy(),reward.copy()])
    states_action_lists.append(states_action_list.copy())
    window_state_lists.append(window_state_list.copy())
    drone_state_lists.append(drone_state_list.copy())
    states_action_list.clear()
    window_state_list.clear()
    drone_state_list.clear()

    # Mean of last 40 episodes
    avg_reward = np.mean(ep_reward_list[-40:])
    print("Episode * {} * Avg Reward is ==> {}".format(ep, avg_reward))
    print("Episode * {} * Reward is ==> {}".format(ep, episodic_reward))
        
    avg_reward_list.append(avg_reward)
    if en_noise_decay:
        ou_noise.std_dev[0]=init_std_dev * math.exp(-(ep/total_episodes)*final_exponent_val)


import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline

# Plotting states
index_of_episode = 0
pos_diff_norm = []
vel_diff_norm = []
acc_diff_norm = []
act_norm = []
act_vals_x, act_vals_y, act_vals_z = [], [], []
inst_reward, accum_reward = [], []
time_vals = np.linspace(0, EnvUAV_cmd_z.episode_timeout, len(states_action_lists[index_of_episode]))

for tuble in states_action_lists[index_of_episode]:
    pos_val_x, pos_val_y, pos_val_z = tuble[0][:3]
    pos_diff_norm.append(np.linalg.norm([pos_val_x, pos_val_y, pos_val_z]))
    
    vel_val_x, vel_val_y, vel_val_z = tuble[0][3:6]
    vel_diff_norm.append(np.linalg.norm([vel_val_x, vel_val_y, vel_val_z]))
    
    acc_val_x, acc_val_y, acc_val_z = tuble[0][6:9]
    acc_diff_norm.append(np.linalg.norm([acc_val_x, acc_val_y, acc_val_z]))
    
    act_vals = tuble[1]
    act_x = act_vals[0]
    act_y = act_vals[1]
    act_z = act_vals[2]
    act_vals_x.append(act_vals[0])
    act_vals_y.append(act_vals[1])
    act_vals_z.append(act_vals[2])
    act_norm.append(np.linalg.norm([act_x, act_y, act_z]))

    inst_reward.append(-(pos_diff_norm[-1] + (vel_reward_weight_factor ** 2) * vel_diff_norm[-1]) * UAV_alt_ctrl.sampling_dt * 100 * 0.1)
    accum_reward.append(sum(inst_reward))

window_vals = window_state_lists[index_of_episode]
drone_vals = drone_state_lists[index_of_episode]


plt.close()
figure(figsize=(24, 10), dpi=80)

plt.plot(time_vals,pos_diff_norm)
plt.plot(time_vals,vel_diff_norm)
plt.plot(time_vals,acc_diff_norm)
#plt.plot(time_vals,act_norm)
#plt.plot(time_vals,act_vals_z)
plt.plot(time_vals,act_vals_x)
plt.plot(time_vals,act_vals_y)
plt.plot(time_vals,act_vals_z)
#plt.plot(time_vals,accum_reward)
#plt.plot(time_vals,inst_reward)


plt.xlabel("Time")
plt.ylabel("State and Action Trajectories")
plt.rcParams.update({'font.size': 20})

plt.legend(['Norm Position','Norm Velocity','Norm Acceleration','Norm Actor','Actor_x','Actor_y','Actor_z','accumelated reward (x0.1)','Instantaneous reward'])
plt.grid()
plt.show()

drone_x = [item[0] for item in drone_vals]
drone_y = [item[1] for item in drone_vals]
drone_z = [item[2] for item in drone_vals]
window_x = [item[0] for item in window_vals]
window_y = [item[1] for item in window_vals]
window_z = [item[2] for item in window_vals]
# Create a figure and axis for the plot
fig, ax = plt.subplots()
plt.plot(time_vals,drone_y)
plt.plot(time_vals,window_y)
ax.set_xlim(0, max(max(time_vals), max(time_vals)))
ax.set_ylim(min(min(drone_y), min(window_y)) - 1, max(max(drone_y), max(window_y)) + 1)
plt.xlabel('Time')
#plt.ylabel('y')
#plt.title('Drone Altitude and Window Altitude vs. Time')
plt.legend(['Drone ','Window '])
plt.grid()
plt.show()
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter3D(drone_x, drone_y, drone_z,'blue');
#ax.scatter3D(window_x, window_y, window_z, 'gray') 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3D
import numpy as np
import time
# Define the YourClassNameHere class
np.savetxt('drone_x.csv', [drone_x], delimiter=',')
np.savetxt('drone_y.csv', [drone_y], delimiter=',')
np.savetxt('drone_z.csv', [drone_z], delimiter=',')
print(drone_x)
class YourClassNameHere:
    def __init__(self,drone_x,drone_y,drone_z,window_x,window_y,window_z):
        self.fig = plt.figure(figsize=(12, 8))
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.counter = 0
        self.drone_x = drone_x
        self.drone_y = drone_y
        self.drone_z = drone_z
        self.window_x = window_x
        self.window_y = window_y
        self.window_z = window_z

    def animate(self, num_frames):
        for frame in range(num_frames):
            self.counter += 1
            #self.ax.cla()

            # Plot the window's 3D position
            window_position = self.ax.scatter(self.window_x[5*frame], self.window_y[5*frame], self.window_z[5*frame], c='r', marker='o')
            drone_position = self.ax.scatter(self.drone_x[5*frame], self.drone_y[5*frame], self.drone_z[5*frame], c='r', marker='o')

            # Set labels and limits for the 3D plot
            self.ax.set_xlabel('X (m)')
            self.ax.set_ylabel('Y (m)')
            self.ax.set_zlabel('Z (m)')
            self.ax.set_xlim([-1, 1])
            self.ax.set_ylim([-1, 1])
            self.ax.set_zlim([-1, 1])
            self.ax.legend()

            # Redraw the 3D plot
            plt.draw()
            plt.pause(0.001)

    def render(self, mode='human', close=False):
        pass  # Leave this empty for animation


animation_example = YourClassNameHere(drone_x,drone_y,drone_z,window_x,window_y,window_z)
# animation_example.animate(num_frames=len(time_vals))
animation_example.animate(len(time_vals))
plt.show()

#print(window_state_lists[index_of_episode])

                    