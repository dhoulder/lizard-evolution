######### enviornmentalChange #########

# Function changes the values of the background environment based on:
# 1 - the model of environmental change (linear or sine wave)
# 2 - the parameters of those models (freq and amp for sine model, slope for linear model)
# 3 - the heterogeneity of environmental change magnitude across domain (binary variable: heterogenous or homogenous)

############# Arguments ###############

# env = current environmental layer
# start.env = original environmental layer (at timestep = 1)
# time = current timestep
# freq = frequency of sine wave (for model = sine)
# amp = amplitude of sine wave (for model = sine)
# slope = slope of linear increase (for model = linear)
# hetero = logical. specifies whether enviro change is spatially homogenous or heterogenous
# env.change.matrix = matrix that matches env domain that specifies how much the environment changes in each grid cell

enviroChange <- function (start.env, env, current.time, amp, freq, slope, model, hetero=T, env.change.matrix) {

# if homogenously varying across domain   
if(hetero==F){
# for linear model environment 
  if(model == "linear") {
    env <-  env + slope 
    } else {
# for sine model environment 
    env <-  start.env + (sin((current.time-1)/freq)/amp) 
    }
} else {
# if hetergoensouly changing across across domain   
  if(model == "linear"){
# for linear model 
    data <- matrix(env@data@values, ncol=100, byrow=F)
    data <- data + (slope*env.change.matrix)
    env@data@values <- as.numeric(data)
  } else {
# for sine model
    data<-matrix(start.env@data@values, ncol=100, byrow=F)
    data <- data + ((sin((current.time-1)/freq)/amp)*env.change.matrix)  
    env@data@values <- as.numeric(data)
  }
} 
#returns new env layer    
return(env)
}




