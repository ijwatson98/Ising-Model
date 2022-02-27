# Ising-Model

To run in terminal:
python 'checkpoint 1_faster.py'
You’ll then be prompted for the number of steps, the dynamic type (Glauber or Kawasaki), the lattice size, the initial temperature, the final temperature, the temperature step, and whether you wish to animate. The submitted graphs are using the following inputs:
python 'checkpoint 1_faster.py' 
10000
G
50
1
3
0.1
N
Code will display animation every 100 steps (if animation selected) and produce graphs of average magnetisation, susceptibility, average energy, and heat capacity for Glauber dynamics. For Kawasaki dynamics it will produce graphs of average energy and heat capacity (magnetisation remains at 0 due to the nature of swaps).
