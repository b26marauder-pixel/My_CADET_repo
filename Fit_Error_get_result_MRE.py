# %%
# Component System

from CADETProcess.processModel import ComponentSystem 
component_system = ComponentSystem()

component_system.add_component ('H2O', charge=0)
component_system.add_component ('H3O+', charge=1)
component_system.add_component ('OH-', charge=-1)

component_system.add_component('IA', species= ['H2IA','HIA-','IA2-'])

# %%
from CADETProcess.processModel import MassActionLaw

k_Hin_i=1.0

reaction_system=MassActionLaw(component_system)
#Wasser Eigendissoziation
reaction_system.add_reaction(
    components=['H2O','H3O+','OH-'],
    coefficients=[-2,1,1],
    k_fwd=k_Hin_i  ,    #
    k_bwd=3.0629747610879686e+17 ,      #
)

#Itaconsäure Dissoziation Stufe 1
reaction_system.add_reaction(
    components=['H2O','H2IA','HIA-','H3O+'],
    coefficients=[-1,-1,1,1],
    k_fwd=k_Hin_i  ,    #
    k_bwd=220328.43247112754 ,      #
)

#Itaconsäure Dissoziation Stufe 2
reaction_system.add_reaction(
    components=['H2O','HIA-','IA2-','H3O+'],
    coefficients=[-1,-1,1,1],
    k_fwd=k_Hin_i  ,    #
    k_bwd=1390178.4266546632 ,      #
)


# %%
# Inlet und Outlet


from CADETProcess.processModel import Inlet, Outlet

inlet = Inlet(component_system, name='inlet')
inlet.flow_rate = 5*1e-6/60  # m^3 / s

inlet.c = [[5.54227005e+04],[9.98705947e-06],[1.00414010e-03],[3.91321182e-07],[9.85627606e-03],[3.93452548e+01]]
 

outlet = Outlet(component_system, name='outlet')

# %%
from CADETProcess.processModel import BiLangmuir

binding_model_Bi = BiLangmuir(component_system)
binding_model_Bi.n_binding_sites=3
binding_model_Bi.is_kinetic = True

binding_model_Bi.adsorption_rate=[0, 0, 0, 1000000.0, 0, 0, 0, 0, 0, 0, 1000000.0, 0, 0, 0, 0, 0, 0, 1000000.0]
binding_model_Bi.desorption_rate=[0, 0, 0, 19810296.599760693, 0, 0, 0, 0, 0, 0, 4846402.951265541, 0, 0, 0, 0, 0, 0, 549028220.0505106]
binding_model_Bi.capacity=[1, 1, 1, 2836.2951575710995, 1, 1, 1, 1, 1, 1, 898.739431206764, 1, 1, 1, 1, 1, 1, 655.6802459646426]

# %%
# Trennsäule

from CADETProcess.processModel import LumpedRateModelWithPores

column = LumpedRateModelWithPores(component_system, name='column')
column.check_required_parameters
column.binding_model = binding_model_Bi
column.bulk_reaction_model=reaction_system

column.length = 0.246 # m
column.diameter = 0.016 # m
column.bed_porosity = 0.3863398
column.particle_porosity = 0.09563789
column.particle_radius = 286.3319*1e-6/2 # m

column.axial_dispersion = 5.81359520865743e-07   # m^2 / s
column.film_diffusion = [1e-10, 5.8e-6, 5.8e-6, 4.5e-6, 5.3e-6, 1.9e-6]

column.c = [55344.0, 0.00010000009999999998, 0.0011, 1e-10, 1e-10, 1e-10]
column.cp = [55344.0, 0.00010000009999999998, 0.0011, 1e-10, 1e-10, 1e-10]
column.q = [1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10, 1e-10]

column.discretization.ncol = 55



# %%
#Flow Sheet

from CADETProcess.processModel import FlowSheet

flow_sheet = FlowSheet(component_system, name='flow_sheet')

flow_sheet.add_unit(inlet)
flow_sheet.add_unit(outlet, product_outlet=True)
flow_sheet.add_unit(column)

flow_sheet.add_connection(inlet, column)
flow_sheet.add_connection(column,outlet)

# %%
#Prozessparameter

from CADETProcess.processModel import Process

process = Process(flow_sheet,'Chromatographie')
process.cycle_time = 160*60

# %%
from CADETProcess.simulator import Cadet

process_simulator = Cadet() 

process_simulator.time_integrator_parameters.abstol=1e-12
process_simulator.time_integrator_parameters.reltol=1e-6

# %%
simulation_results = process_simulator.simulate(process)

# %%
#Laborversuche laden

import pandas as pd
data = pd.read_excel("IA-values_pH8_MRE.xlsx", index_col=0)

data.plot()

# %%
#Reference

from CADETProcess.reference import ReferenceIO

reference = ReferenceIO('experimental data',data.index*60,data,component_system=component_system)

# %%
from CADETProcess.comparison import Comparator
comparator = Comparator()

comparator.add_reference(reference)

#Difference metrics

comparator.add_difference_metric('RMSE', reference, 'column.outlet',components='IA',use_total_concentration_components=True)

#Plot comparison und Differenz ausgeben

comparator.plot_comparison(simulation_results)

#Difference metrics

metrics = comparator.evaluate(simulation_results)
print('Metrics Auswertung:',metrics)

# %%
from CADETProcess.optimization import OptimizationProblem
optimization_problem = OptimizationProblem('film_diffusion')

optimization_problem.add_evaluation_object(process)

#Optimization Problem

###Variablen###

optimization_problem.add_variable(
name='film_diffusion_H2O', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[0]
)

optimization_problem.add_variable(
name='film_diffusion_H3Oplus', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[1]
)

optimization_problem.add_variable(
name='film_diffusion_OHminus', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[2]
)

optimization_problem.add_variable(
name='film_diffusion_H2IA', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[3]
)

optimization_problem.add_variable(
name='film_diffusion_HIAminus', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[4]
)

optimization_problem.add_variable(
name='film_diffusion_IA2minus', parameter_path='flow_sheet.column.film_diffusion',
lb=1e-9, ub=1e-3,
transform='auto',
indices=[5]
)

#Evaluator und Objectives

optimization_problem.add_evaluator(process_simulator)

optimization_problem.add_objective(
comparator,
n_objectives=comparator.n_metrics,
requires=[process_simulator]
)

#Optimizer-Einstellungen

from CADETProcess.optimization import U_NSGA3
optimizer = U_NSGA3()

optimizer.n_cores = 8   
optimizer.pop_size = 32
optimizer.n_max_gen = 8

def callback_funcion(simulation_results, individual):
    print('Individual:', individual)

optimization_problem.add_callback(callback_funcion, requires=[process_simulator])



# %%
#Optimieren

optimization_results = optimizer.optimize(
    optimization_problem,
    use_checkpoint=False
    )

# %%
print(optimization_results.x)
print(optimization_results.f)

optimization_results.plot_convergence('objectives')

optimization_results.plot_objectives(autoscale=False)


