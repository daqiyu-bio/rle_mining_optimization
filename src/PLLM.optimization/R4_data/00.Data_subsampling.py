import pandas as pd
import numpy as np

data = pd.read_csv( "01.R4_A43_SingleMutation_chosen_modeling_20250705.csv" )
data = data[ [ "Order", "Mutation", "Mean", "input", "Class" ] ]
data.columns = [ 'Order', 'input_ids', 'labels', 'input', 'Class' ]

mean_value     = np.mean(data.query( 'Class == "Train"' )['labels'])
std_value      = np.std(data.query( 'Class == "Train"' )['labels'])
data['labels'] = data.apply( lambda x:(x['labels'] - mean_value) / std_value, axis = 1 )

for seed in [101, 102, 103, 104, 105]:
    data2 = data.query( 'Class == "Train"' ).copy()
    data_subsample1 = data2.sample( n = 80, random_state = seed ).copy()
    data2[ 'Class' ] = data2.apply( lambda x: "Train" if x['input_ids'] in list(data_subsample1['input_ids']) else "Valid", axis = 1 )
    data2.to_csv( ".".join( [ "01.R4_A43_SingleMutation_chosen_modeling_20250705.subsampled", str(seed), "csv" ] ) )

