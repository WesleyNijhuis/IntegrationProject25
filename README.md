# IntegrationProject25
In this Integration project, we aim to design multiple controllers for a furuta pendulum. A report of the project can be found in "Identification and control of a furuta pendulum.pdf".

The project consisted of three fases. First, a first principles model of the Furuta pendulum was established and equipment was callibrated. Next, the model parameters were identified through identification around the stable equilibirum. Uniquely to the approach we used is the fact that subspace identification methods were used in combination with a structuring state variable transformation and a refining OE identification cycle(more information in chapter 4 and the appendix). In the third phase, a Kalman observer, LQR, swing-up controller and model predictive controller were synthesized and implemented. For the MPC, an optimized matlab s-function was created to achieve the speed necessary for a long prediction horizon. 


## The Furuta pendulum
![pendulum](./Report/figures/furutapendulum.png)

## Schemes of LQG and MPC
![mpc](./Report/figures/MPC.png)

![mpc](./Report/figures/LQG.png)

## Results of OE enhancement vs PEM enhancement after subspace identification and strucuring transformation
![mpc](./Report/figures/OEvsPEM.png)

## Terminal set of MPC controller
![Xf](./Report/figures/X_f.png)

## Tracking results of MPC in the presence of disturbances
![tracking](./Report/figures/Demo_MPC.png)





