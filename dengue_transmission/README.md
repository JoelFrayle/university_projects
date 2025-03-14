# Mathematical Model for Dengue Transmission Dynamics with Seasonal Variations and Control Strategies  

This repository contains a project from my undergraduate work in Biomedical Engineering, where we developed a mathematical model to simulate the transmission dynamics of dengue virus in human and mosquito populations. The model is based on the epidemiological context of Medellín, Colombia, and evaluates the effectiveness of various control strategies (and their combinations) to determine the most effective approach for reducing dengue transmission.  

## Overview  
- **Paper**: The [`Dengue_Transmission_Model.pdf`](./Dengue_Transmission_Model.pdf) file provides a detailed academic paper (originally written in Spanish, now translated to English) with:  
  - Biological background on dengue virus, its transmission, and the role of *Aedes aegypti* mosquitoes.  
  - Mathematical modeling approach (system of differential equations, parameters, and seasonal variations).  
  - Simulation results, comparison of control strategies, and potential future work.  

- **Python Script**: The [`dengue_transmission_model.py`](./dengue_transmission_model.py) implements the mathematical model. It simulates:  
  - Human populations (susceptible, exposed, infected, and recovered).  
  - Mosquito populations (susceptible, exposed, and infected).  
  - Seasonal variations (rainy, normal, and dry seasons).  
  - Control strategies (insecticide use, larvicide application, mechanical control, and vaccination).  

## Key Features  
- **Differential Equations**: Continuous modeling of population changes in humans and mosquitoes, incorporating seasonal variations.  
- **Control Strategies**: Evaluation of individual and combined control methods to determine their effectiveness in reducing dengue transmission.  
- **Visualization**: Graphs showing the impact of control strategies on human and mosquito populations over time.  
- **Seasonal Variations**: Incorporation of seasonal changes to reflect real-world conditions in Medellín, Colombia.  

## How to Use  
1. Clone this repository:  
   ```bash  
   git clone https://github.com/your-username/dengue-transmission-dynamics.git  
