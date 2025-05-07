# Benchmark outputs for Q8

These files were created with runbenchq8.sh
OMP_NUM_THREADS unset.
They contain the benchmark table and first 50 lines of profiler output.
Profiler stats are sorted by CUMULATIVE time. Files with ".time" extension are sorted by TOTAL time.

## Profiler results for GBTFITSLOAD
These are the top CPU using methods as listed by the profiler

### with flags
  
   - Selection._check_for_duplicates
     - this method checks that a selection rule (DataFrame) does not create an indentical DataFrame to an existing selection rule.   It uses DataFrame.equals().   Because every flag rule is stored as a DataFrame, this is called N-1 times for the Nth flag rule.   There may be savings available by finding abetter way to compare DataFrames (e.g. assert_frame_equals which compares column by column or some more clever way)
 - Selection._addrow
    - This calls astropy Table.add_row which can be expensive. There may be a better way to build the "flag/selection summary table" other than row by row when we know we will be reading in multiple rows of flags.
- Selection._add_utc_column 
    - We add the UTC column so that users can easily select data by time using a Time object. 
    - Despite using the "array version" of astropy Time instantiation, this call is still expensive when the number of rows is large.   Pandas supported np.datetime64, which may be a better way to store this column.   And just do a conversion of the user's desired time range to datetime64
   - Time spent in millions of calls to astropy Table functions
 

### without flags

- Selection._add_utc_column 
     - see above. Note for some reason this method is called 8 time regardless of how many files are loaded.  (Same for with flags)

- Various astropy and pandas methods. Some astropy methods appear to be related to initial manipulation of the binary table (drop_fields) when creating the index.  
 
