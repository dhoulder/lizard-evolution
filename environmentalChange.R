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

  if (amp==0) {
    return(start.env)

  } else {

    # if homogenously varying across domain
    if(hetero==F){
    # for linear model environment
      if(model == "linear") {
        env <-  env + slope
        } else {
    # for sine model environment
        env <-  start.env + (sin((2 * pi * (current.time-1)) / freq) * amp)
        }
    } else {
    # if hetergenously changing across across domain
      if(model == "linear"){
    # for linear model
        data <- matrix(env@data@values, ncol=100, byrow=F)
        data <- data + (slope*env.change.matrix)
        env@data@values <- as.numeric(data)
      } else {
    # for sine model
        env <- start.env # this creates a raster of the same dimensions - all the values will be replaced
        env[] <- start.env[] + ((sin((2 * pi * (current.time-1)) / freq) * amp) * env.change.matrix[])
      }
    }

  #return new env layer
  return(env)
  }
}




