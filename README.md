---
typora-copy-images-to: upload
---

# Computer Animation and Special Effects HW1

## Build on Microsoft Windows with Visual Studio 2017/2019

### Instruction

- Open SoftSim.sln
- Build
- Executable will be in ./bin

## Build on other platforms and/or compilers

### :warning: **This method is not well-tested, so it may not work properly.**

### :warning: **Thus, you are expected to solve buggy or trivial problems yourself.**

### Some tested platforms (cmake 3.19):

- Ubuntu Groovy Gorilla (20.10) with GCC   10.2.0
- macOS Catalina        (10.15) with Clang 11.0.3

### Prerequisite

- [Cmake](https://cmake.org) (version >= 3.14)
- Compiler (e.g. GCC)

### Instruction

- Run:
```bash=
cmake -S . -B build
cmake --build build --config Release --target install --parallel 8
```
- Executable will be in ./bin

### If you are building on Linux, you need one of these dependencies, usually `xorg-dev`

- `xorg-dev` (For X11)
- `libwayland-dev wayland-protocols extra-cmake-modules libxkbcommon-dev` (For Wayland)
- `libosmesa6-dev` (For OSMesa)



# Soft Sim 110550047

## Introduction/Motivation

This project practice four integrating methods to simulate particles' behaviors. For the particle object, I construct a complex net struct of soft sim. It contains three connection types of spring, **stretch** (connect to neighbor), **shear** (connect to surrounding) and **bending** (connect to the neighbor of neighbor). In each spring struct, we take two endpoints and treat them as particles, and the inner force is respect to **spring force** and **damper force**. For the integrating method, we implement **Explicit Euler method**, **Implicit Euler method**, **Midpoint Euler method** and **Ronge Kutta Fourth method**. These methods simulate the next motion of each particles, depending on the force field, coefficient, terrain, inner force.

## Fundamentals

### Spring

<img src="https://i.imgur.com/AKX09RM.png" style="zoom:80%;" />

+ Stretch : Blue Line
+ Shear : Yellow Line
+ Bending : Green Line
+ Each type of connection has their own spring coefficient and damper coefficient.

### Force

+ Spring Force
  $\vec{f}_a = -k_s(|\vec{x_a} - \vec{x_b}|-r)\frac{\vec{x_a}-\vec{x_b}}{|\vec{x_a}-\vec{x_b}|}$

+ Damper Force
  $\vec{f}_a = -k_d(|\vec{v_a} - \vec{v_b}|)\cdot(|\vec{x_a} - \vec{x_b}|)\frac{\vec{x_a}-\vec{x_b}}{|\vec{x_a}-\vec{x_b}|^2}$

+ Internal Force
  Update two endpoint and be care of the force direction.

  ```cpp
  particles[spring.getSpringStartID()].addForce(springForce + damperForce);
  particles[spring.getSpringEndID()].addForce(-springForce - damperForce);
  ```

### Terrain

+ Collision Detection
  $$
  ||N|| = 1\\
  \vec{N}\cdot(\vec{x}-\vec{p}) < \epsilon\\
  \vec{N}\cdot\vec{v} < 0
  $$

+ Collision Response : Quick modification
  $$
  v' = -k_rv_N + v_T
  $$

+ Contact Detection
  $$
  ||N|| = 1\\
  |N\cdot(x-p)| < \epsilon\\
  |N\cdot v| < \epsilon\\
  $$

+ Contact Forces : Add Force
  $$
  f^c = -(N\cdot f)N\  \ \ if\ N\cdot f < 0\\
  f^f = -k_f(-N\cdot f)v_t \  \ \ if\ N\cdot f < 0
  $$

### Integrator

+ **Explicit Euler method**

  $x_{n+1}$ = $x_n$ + $\Delta V_n$

  $v_{n+1}$ = $v_n$ + $\Delta A_n$

  ```c++
  jelly->getParticle(i).addPosition(jelly->getParticle(i).getVelocity() * particleSystem.deltaTime);
  jelly->getParticle(i).addVelocity(jelly->getParticle(i).getAcceleration() * particleSystem.deltaTime);
  ```

+ **Implicit Euler method**

  $x_{n+1}$ = $x_n$ + $\Delta V_{n+1}$

  $v_{n+1}$ = $v_n$ + $\Delta A_{n+1}$

  ```cpp
  particleSystem.computeJellyForce(*jelly);
  for (int i = 0; i < jelly->getParticleNum(); i++) {
              jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
  }
  ```

  

+ **Midpoint Euler method**

  $X_{n+\frac{1}{2}}$ = $x_n$ + $\frac{1}{2}\Delta V_{n+1}$

  $V_{n+\frac{1}{2}}$ = $v_n$ + $\frac{1}{2}\Delta A_{n+1}$s

  $x_{n+1}$ = $x_n$ + $\Delta V_{n+\frac{1}{2}}$

  $v_{n+1}$ = $v_n$ + $\Delta A_{n+\frac{1}{2}}$

  ```cpp
  for (int i = 0; i < jelly->getParticleNum(); i++) {
              origin.push_back(jelly->getParticle(i));
              //Eigen::Vector3f deltaX = 0.5 * particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
              Eigen::Vector3f eulerV = jelly->getParticle(i).getAcceleration() * particleSystem.deltaTime;
              Eigen::Vector3f eulerP = jelly->getParticle(i).getVelocity() * particleSystem.deltaTime;
              jelly->getParticle(i).addVelocity(eulerV * 0.5);
              jelly->getParticle(i).addPosition(eulerP * 0.5); 
              jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
  }
          particleSystem.computeJellyForce(*jelly);
  for (int i = 0; i < jelly->getParticleNum(); i++) {
              jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
              jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
  }
  ```

  

+ **Ronge Kutta Fourth method**
  $$
  x(t_0+h) = x(t_0) + \frac{1}{6}(k_1+2k_2+2k_3+k_4)\\
  v(t_0+h) = v(t_0) + \frac{1}{6}(k_1+2k_2+2k_3+k_4)\\
  \\
  k_1 = hf(x_0,t_0)\\
  k_2 = hf(x_0+\frac{k_1}{2}, t_0 + \frac{h}{2})\\
  k_3 = hf(x_0+\frac{k_2}{2}, t_0 + \frac{h}{2})\\
  k_4 = hf(x_0 + k_3,t_0+h)\\
  $$

  ```cpp
  for (int i = 0; i < jelly->getParticleNum(); i++) {
              StateStep k;
              k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
              k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
              jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity() / 2);
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration() / 2);
              k1.push_back(k);
              
          }
          particleSystem.computeJellyForce(*jelly);
          for (int i = 0; i < jelly->getParticleNum(); i++) {
              StateStep k;
              k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
              k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
              jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity() / 2);
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration() / 2);
              k2.push_back(k);
          }
          particleSystem.computeJellyForce(*jelly);
          for (int i = 0; i < jelly->getParticleNum(); i++) {
              StateStep k;
              k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
              k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
              jelly->getParticle(i).setPosition(origin[i].getPosition() + particleSystem.deltaTime * jelly->getParticle(i).getVelocity());
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + particleSystem.deltaTime * jelly->getParticle(i).getAcceleration());
              k3.push_back(k);
          }
          particleSystem.computeJellyForce(*jelly);
          for (int i = 0; i < jelly->getParticleNum(); i++) {
              StateStep k;
              k.deltaVel = particleSystem.deltaTime * jelly->getParticle(i).getAcceleration();
              k.deltaPos = particleSystem.deltaTime * jelly->getParticle(i).getVelocity();
              k4.push_back(k);
          }
          for (int i = 0; i < jelly->getParticleNum(); i++) {
              jelly->getParticle(i).setPosition(origin[i].getPosition() + (k1[i].deltaPos + 2 * k2[i].deltaPos + 2 * k3[i].deltaPos + k4[i].deltaPos) / 6);
              jelly->getParticle(i).setVelocity(origin[i].getVelocity() + (k1[i].deltaVel + 2 * k2[i].deltaVel + 2 * k3[i].deltaVel + k4[i].deltaVel) / 6);
              jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
          }
  ```

  

## Result and Discussion

### The difference between integrators 

Below Screenshot will keep this format.

| Explicit Euler | Implicit Euler  |
| -------------- | --------------- |
| **Midpoint**   | **Ronge Kutta** |

+ First Bounce

  **Midpoint** and **Ronge Kutta Fourth** method show the more deformation value at local range, while **Explicit Euler** and **Implicit Euler** method look like alleviating the impact force with the global range.

  <img src="https://i.imgur.com/TnPERwS.png" />![](https://i.imgur.com/k5XQiop.png) 
  ![](https://i.imgur.com/MNe2Laz.png)![](https://i.imgur.com/Rc4Oimi.png)

+ Slide Angle
  Because the previous difference, the angle of sim to slide is distinct. **Explicit Euler** and **Implicit Euler** slide with its edge, while **Midpoint** and **Ronge Kutta Fourth** slide with its face.

  ![](https://i.imgur.com/3gmnYaV.png)![](https://i.imgur.com/YldJuzd.png)

![](https://i.imgur.com/4W0QQGw.png)![](https://i.imgur.com/cwwNQIm.png)

+ First Touch Margin

  **Explicit Euler** and **Implicit Euler** slide out and quickly fall back, the **Midpoint** one stay a bit longer time then falling back, and the **Ronge Kutta Fourth** almost stay in there but eventually fall back to hole.
  ![image-20230325113130169](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113130169.png)![image-20230325113150425](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113150425.png)
  ![image-20230325113227799](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113227799.png)![image-20230325113106514](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113106514.png)

+ Second Touch Margin

  The previous one apply different force to make sim rotate around y-axis. The **Explicit Euler** one seems stably, while the **Implicit Euler** one turns more angle and suffer from more contact angel so that hardly touch the same height as previous one.

  The **Midpoint** turns approximately 30 degree, and the **Ronge Kutta Fourth** turns almost 90 degree.

  ![image-20230325113450928](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113450928.png) ![image-20230325113520815](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113520815.png)
  ![image-20230325113558937](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113558937.png) ![image-20230325113635958](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325113635958.png)

### Effect of parameters 

+ Spring Coefficient 
  1000、2000、4000、8000, obviously bounce higher
  ![image-20230325120119768](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120119768.png) ![image-20230325115929486](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325115929486.png) 
  ![image-20230325115949268](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325115949268.png) ![image-20230325120007706](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120007706.png)

+ Damper Coefficient
  10

  ![image-20230325120526028](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120526028.png) ![image-20230325120559727](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120559727.png)
  60

  ![image-20230325120348301](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120348301.png) ![image-20230325120634715](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325120634715.png)

+ Resist Coefficient
  0.8 : ![image-20230325134249227](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134249227.png)0.2 : ![image-20230325134155176](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134155176.png)

+ Friction Coefficient

  0.2 : ![image-20230325134249227](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134249227.png)0.8 : ![image-20230325134334649](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134334649.png)

+ Jelly Rotation

  Deformation Duration
  <img src="C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134430173.png" alt="image-20230325134430173" style="zoom: 80%;" /> ![image-20230325134453798](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134453798.png)

  ![image-20230325134513810](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134513810.png)![image-20230325134543625](C:\Users\wujoe\AppData\Roaming\Typora\typora-user-images\image-20230325134543625.png)

## Conclusion

In this project, we simulate many force and influence factors to show the difference between integrate methods. In general, adjust the spring coefficient and damper coefficient much lower and enhance friction coefficient will simulate just like a real sim there. Depending on situations, we can choose different method to  implement.
**Explicit** is much simple but not seems real at local part.
**Implicit** also not seems real at local part, but is much reasonable in this project.
**Midpoint** can observe the deformation and is the most suitable method in this method.
**Ronge Kutta** can observe the deformation at local part, but it may be much stubborn and heavy when simulating the sim.
Above all, not every situation has the commonly best answer. If so, we also need to simulate every part as a real world, such as the small force inside, but it is not essential and necessary. Under the simulation, we can trial and error and see the much closer to the reality or imagine one.
