# Scheduler for the Xamidimura telescopes
An heuristic dispatch scheduler which selects the current best eclipsing binary target to observe.

## Installation
Make sure the entire directory is downloaded before attempting to run the scheduler.

### Non-standard packages you will need
- astroplan
- prettytable
- sqlite3

## Usage

There are three schedulers included in the package: 
1. `scheduler.py` --- User-friendly interface, most function-rich of the three
2. `dumb_scheduler.py` --- Most simple scheduler, see below for example
3. `weather_sim_scheduler.py` --- Specifically for running simulations on weather logs

### 1. `scheduler.py`

To run the scheduler with an interface, in scheduler.py ensure 
```python
interface = True
```
then run the scheduler script in terminal.

To run the scheduler without the interface, in scheduler.py ensure
```python
interface = False
```
then use the functions, for example:
```python
import scheduler

# Best target to observe now
best_target, status = scheduler.observe_now()
target_name, phase = best_target.name, best_target.phase

# All targets currently feasible
feasible_targets = scheduler.feasible_now()

# Detailed information about an eclipsing target
scheduler.eclipse_info(target_name, display=True)

# Information on why a specific target is not feasible
scheduler.why_not(target_name)

# Tonight's best targets
target_names, target_info = scheduler.tonights_best()

# Sky plot of target(s) over the night's duration
scheduler.make_sky_plot(target = best_target)
scheduler.make_sky_plot(target_array = target_names)
```

### 2. `dumb_scheduler.py`

Run the script and use the functions as follows:
```python
# Make target list from database - only need to do once
target_list = ds.make_eb_objects(dbconn, moon_tol)

# Initially not observing a target
current_target = None

# Selects best target at current time and retrieve info - can loop this
time = Time.now()
current_target = observe_now(time, target_list, current_target)
target_name, phase = current_target.name, current_target.phase
```

### 3. `weather_sim_scheduler.py`

Just run in terminal and wait (a progress bar will appear). Outputs result to text files in `weather_sim/` subdirectory.
Can select observing intervals of either an hour or 15 minutes. In `main()`:
```python
mode = 'hourly' #'quarterly'
```