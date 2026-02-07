Full code will be available after publication of the paper

Files that will be added are:

- plot_settings.m
- plot_functions.m
  
These files can be used to get the outputs similar to the paper.

# Stimulation of a Neuron-Astrocyte Network (SNAN) with Working Memory (WM)

Some parts of the code are from: https://github.com/PedRaMNG/working-memory-with-self-repairing

This repository contains MATLAB code for online/offline excitatory/inhibitory stimulation of a Neuron-Astrocyte Network (SNAN) with Working Memory (WM), published in iScience journal.

Link to the paper: ...

# Requirements

* The code was tested on MATLAB 2022b, it might work with earlier versions.
* The minimum required amount of RAM is 32 GB.

# Settings

## Quick setup

There are 3 experiments in the paper, each requiring a unique set of settings. To simplify this for users, we have grouped the variables that need to be changed in a switch statement inside `model_parameters.m`. All you need to do is change the `params.simPattern` variable at the top of the `model_parameters.m` file.

```matlab
params.simPattern = 1; % Experiment 1 - Enhancement (default: damaged network)
params.simPattern = 2; % Experiment 2 - Disruption (default: disrupts both train and test)
params.simPattern = 3; % Experiment 3 - Augmentation
```

* To start the simulation, run the `main.m` file. When you run this file, you will get a GUI pop-up to confirm your settings before start of the simulation. You have to press **continue** to start the simulation otherwise the simulation **will not start**.
* To plot the outputs, open the `plot_settings.m` file, set the variables you want to plot to **1**, and run the file.
* For video, when you get the result, in the video window, go to "playback" -> "frame rate" and increase the frame rate to `500` or `1000` and check the "allow frame drop ..." box. This increases the speed of video playback.
* In the `plot_settings`, make sure that the variable you are trying to plot is with the dimension of params.n in the `init_model.m`.

## Experiment 1 - Enhancement

To get the results for this experiment, set the following parameters in the **`model_parameters.m`** file:

```matlab
params.simPattern = 3;
params.impairmode = 1; % For healthy network.
params.impairmode = 3; % For damaged network (default).
```

**Configuration Details:**

Setting `simPattern` to **3** activates the **Enhancement** mode (Case 3). This configuration:

* Sets the simulation duration (`t_end`) to **5.6s**.
* Enables **Synaptic Damage Mode 3** (`impairmode = 3`), which applies concentrated damage based on a specific damaging pattern.
* Sets the learning and testing impulse order to `[1,2,3,4]`.
* Enables **Automatic Stimulation** (`stimulation_mode = 1`) targeting **Excitatory neurons** (`stimulation_neuro_E = 1`) at **3.2s** to counteract any damage. Activating these two parameters uses the `stimulator_enhance` function in `simulate_model.m`.

Since this experiment involves two parts, a healthy network and a damaged network, you can enable synaptic damage by setting `impairmode = 3`, which applies concentrated damage based on a specific damaging pattern. For a healthy network, set `impairmode = 0`. In healthy mode, even though automatic stimulation is active, the Enhancement algorithm recognizes this and does not activate any stimulating electrodes.

## Experiment 2 - Disruption

To get the results for this experiment, set the following parameters in the **`model_parameters.m`** file:

```matlab
params.simPattern = 4;
params.disruptBothTrainTest = 0; % Disrupts only in the training phase.
params.disruptBothTrainTest = 1; % Disrupts both in train and test phases (default).
```

**Configuration Details:**
Setting `simPattern` to 4 activates the **Disruption** mode (Case 4). This configuration:

* Sets the simulation duration (`t_end`) to **3.5s**.
* Sets the learning and testing impulse order to `[9, 3, 11, 7]`.
* Maintains a **Healthy Network** (`impairmode = 0`) but applies external disruption.
* Sets disruption pattern indices to `[3,7]`. To disrupt other patterns, change the `disruption_pattern` variable.
* Enables **Automatic Stimulation** (`stimulation_mode = 1`) targeting **Inhibitory neurons** (`stimulation_neuro_I = 1`) at **1.28s** to disrupt the pattern recall. Activating these two parameters uses the `stimulator_dist` function in `simulate_model.m`.

## Experiment 3 - Augmentation

To get the results for this experiment, set the following parameters in the **`model_parameters.m`** file:

```matlab
params.simPattern = 5;
```

**Configuration Details:**
Setting `simPattern` to **5** activates the **Augmentation** mode (Case 5). This configuration:

* Sets the simulation duration (`t_end`) to **2.6s**.
* Sets learning order to `[4, 5]` and testing order to `[4, 2, 5, 1]`.
* Maintains a **Healthy Network** (`impairmode = 0`).
* Switches to **Manual Stimulation** (`stimulation_mode = 0`) with **Augmentation enabled** (`stimulation_Augmentation = 1`).
* Activates specific electrodes `[15,8,12,16]` for stimulation.
* Uses two sets of timing for manual stimulation (`stimulation_time_neuro_E_Augmentation = [0.12, 0.42]`).
* Activates specific electrodes for each timing `[6,10,14]` and `[3,4,11,12]`.
* Unlike other modes, there is no separate stimulator function for this mode. Instead, the setup is configured inside the `init_model.m` file.

---

# Code

* The term `train` code is mentioned as **sample phase** in the paper.
* In all modes, stimulation current (`I_stim`) is added to `I_app_astro` inside the `get_neuron_astrozone_activity.m` function.
* Files with the **TB** extension are **Test Bench** files used for fast iteration and testing of a function. They are typically designed for standalone execution but sometimes they need data from a simulation. To get the UT files for enhance and disrupt functions, email the corresponding author (Mahmood Amiri).
* `accuracy_elec` measures the correlation, mentioned in the paper.
* Extra code for synaptic plasticity is provided in `step_neurons.m` for future work.
* To load your saved data, open `saveLoadData.m` and run the **Import** section.
* To use `sendTelegramMessage` function you need to set some configuration inside its function.
* To understand the relationship between code blocks, we recommend to check out the AI-generated overview by **DeepWiki**, linked at the bottom of this readme file. Since it is a fully automated process by AI, **it can make mistakes**.

# Authors

* Pedram Naghieh - *Software, Formal analysis, Writing – original draft, Methodology* - [PedRaMNG](https://github.com/PedRaMNG)
* Mahmood Amiri - *Conceptualization, Methodology, Supervision, Writing review & editing*
* Herbert Peremans - *Conceptualization, Writing review & editing*

# Cite

P. Naghieh, M. Amiri, and H. Peremans, “**Targeted Neuronal Stimulation Modulates Working Memory via Astrocytic Interactions: Insights from a Spiking Neuron-Astrocyte Model**”, iScience, p. ... , jan. 2025, doi: ...

# License

This project is licensed under the MIT License - see the LICENSE.md file for details.

# DeepWiki

AI-generated overview of the code and the relationship between the code blocks:
