# Scheduler for the Xamidimura telescopes
An heuristic dispatch scheduler which selects the current best eclipsing binary target to observe.

## Installation
Make sure the entire directory is downloaded before attempting to run the scheduler.

### Non-standard packages you will need
- astroplan
- sqlite3
- prettytable

## Usage

There are three schedulers included in the package, plus a handy file for database manipulation: 
1. `scheduler.py` --- Most simple scheduler
2. `scheduler_user.py` --- User-friendly interface, most function-rich of the three
3. `scheduler_weather_sims.py` --- Specifically for running simulations on weather logs
4. `easy_database.py` --- Easily modify contents of target info and priority tables

### 1. `scheduler.py`

Returns the name of the best eclipsing binary target to observe at the current time.
Simply run the script; press `ENTER` to call the `observe_now()` function as many times as you like, or `q` to terminate the script.

### 2. `scheduler_user.py`

This is the most comprehensive of the schedulers. Instead of just blindly requesting a target, you can select what information you want to see and how you want it presented. It can be run with or without the interface.

To run the scheduler with an interface, in `scheduler_user.py` ensure 
```python
interface = True
```
then simply run the script.

To run the scheduler without the interface, in `scheduler_user.py` ensure
```python
interface = False
```
then you can use the functions, for example:
```python
import scheduler_user as s

# Best target to observe now
best_target, status = s.observe_now()
target_name, phase = best_target.name, best_target.phase

# All targets currently feasible
feasible_targets = s.feasible_now()

# Detailed information about an eclipsing target
s.eclipse_info(target_name, display=True)

# Information on why a specific target is not feasible
s.why_not(target_name)

# Tonight's best targets
target_names, target_info = s.tonights_best()

# Sky plot of target(s) over the night's duration
s.make_sky_plot(target = best_target)
s.make_sky_plot(target_array = target_names)
```

### 3. `scheduler_weather_sims.py`

Adjust preferences (e.g. score weightings in `EclipsingBinary.mean_priority()` or parameters at the start of `main()`), then run the file and wait --- a progress bar will appear. The simulation output will be saved to text files in the `weather_sim/` subdirectory.

The contents of `feasibility_and_priority.py` have been reproduced in `scheduler_weather_sims.py`, which means you can modify anything you like without affecting the regular schedulers.

You can select observing intervals of either an hour or 15 minutes. In `main()`:
```python
mode = 'hourly' #'quarterly'
```
Be aware that changing mode to quarterly can increase the run time from approx. 5h30 to 27 hours.

### 4. `easy_database.py`

Run file, then use the following lines:
```python
dbconn, dbcurs = connect_database.connect_to()
t = targets_to_dataframe(dbconn)
s = priorities_to_dataframe(dbconn)
```
Before adding any new targets, make sure that they are written in the exact same format as `master_targets.csv`.
