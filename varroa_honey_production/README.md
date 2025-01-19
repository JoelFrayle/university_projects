# Mathematical Model for Honey Production Variation in a Honeybee Colony Infected by *Varroa destructor*

This repository contains a project from my undergraduate work in Biomedical Engineering, where we developed a **mathematical model** to assess how *Varroa destructor* infestation impacts the honey production in a honeybee colony (*Apis Mellifera*). The model compares a healthy colony to an infected one, quantifying the percentage drop in honey yield over a 90-day period.

## Overview

- **Paper**: The [`doc/paper_english.tex`](./doc/paper_english.tex) file provides a detailed academic paper (originally written in Spanish, now translated to English) with:
  - Biological background on honeybees and the *Varroa destructor* parasite.
  - Mathematical modeling approach (differential equations, parameters).
  - Simulation results, discussion, and potential future work.

- **Python Script**: The [`src/model_varroa.py`](./src/model_varroa.py) implements the mathematical model. It simulates:
  1. **Larvae (healthy vs. infected)**  
  2. **Adult bees (healthy vs. infected)**  
  3. **Parasites (female and male)**  
  4. **Honey production**  

  and calculates the percentage difference in honey yield between a healthy colony and one infected by *Varroa*.

### Key Features

- **Differential Equations**: Continuous modeling of population changes in healthy/infected bees, larvae, and parasites.  
- **Honey Production Tracking**: Distinguishes contributions from both healthy and infected bees.  
- **Collapse Point**: Observes how, around a certain day, the infected colony hits a collapse phase leading to a drastic drop in population and honey production.
