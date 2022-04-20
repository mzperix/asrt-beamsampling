import sys
sys.path.append('..')
import filters
import numpy as np

data = dict(rt=np.array([123,101,76,99,45,11,93,140]),
            correct_response=[0,1,1,0,1,0,1,1],
            block = [1,1,1,2,2,2,2,2],
            trial = [1,2,3,1,2,3,4,5],
            )

ini = dict(FILTERS = dict(rt_min=75, 
                          rt_max = 110,
                          skip_start_at_block = 2,
                          )
          )

filtered_data = filters.filter(data, ini)
assert np.prod(filtered_data == np.array([76, 93]))
assert np.prod(filters.new_stream(data) == np.int8([1,0,0,1,0,0,0,0]))
print('Filter tests [ OK ]')