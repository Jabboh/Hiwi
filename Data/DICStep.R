###################WAIC Step Functions: 1. Dropping entire covariates and their
###################qassociated quadratic and interaction terms and 2. 
###################Dropping Interactions and quadratic terms
#x:basic model which is supposed to be reduced
#z: names of the covariates, that will be dropped iteratively
library(mgsub)
x <- m_com
z <- c("Unit", "HIDRO_MEAN_500", "NDVI_2000", "DIST_ARROZ", "DIST_urbano")
step_waic_entCov <- function(form, x, z){
  #define the base model to which the reduced models will be compared
  base <- x
  #Create variable w_com that serves as the condition for the while-loop
  #After the first run, it is an object that stores the WAIC-comparisons
  #(It is WAIC-ranking of the compared models)
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #Create a while loop, that keeps dropping variables until the more complex model
  #performs better than the reduced models
  
  #Define gjam model settings for all the refitting 
  types <- c("CA", "CA", "CA")
  #define model/algorithm parameters: 4000 gibbs steps + burnin of 1000
  ml   <- list(ng = 4000, burnin = 1000, typeNames = types)
  while(rownames(w_com)[1] != "base"){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(z) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- base$fit$DIC
    #Run the models (removing one covariate and all its associations in each run)
    #So, we iterate over all the possible covariates in a for-loop
    for (i in z){
      #index
      j <-match(i,z)
      #We run the model without the terms associated with z (linear, quadratic
      #and interaction terms and store the fitted model in st)
      
      #Delete all the terms associated with covariate i in the model formula
      formul <- mgsub(as.character(form)[2], pattern = c(paste0(i, " \\+ "), paste0(" \\+ ", i), paste0(" \\+ I\\(", i, "\\^2\\)")),
                       replacement = c("","",""))
      #make new xdata for that formula
      model_xdata <- stats::model.matrix(as.formula(paste0("~",formul)), m_com$inputs$xdata)    
      x_train <- as.data.frame(model_xdata)
      #make the new model formula in the form gjam requires
      form_gj <- paste(names(x_train)[-1], collapse ="+")
      form_gj <- paste0("~", form_gj)
      #run the gjam model
      st[[j]] <- gjam(form_gj,
                    ydata = y_train, xdata = x_train, modelList = ml)
      #############BIN HIER!########## 
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    #the model named "IA_500" is the model without all the "IA_500" terms
    n <- c(z, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs, The first row contains the preferred model
    #(the one with the lowest WAIC)
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The variable ", rownames(w_com)[1], " and its associations will be dropped"))
    #deleting the variable that we throw out according to WAIC in the variable
    #names vector z
    z <- z[z != rownames(w_com)[1]]
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[1])]]
  }
  #Return the final model
  return(base)#for some reason this does not work
}

#Implementing Step 2 (Dropping interaction terms) and 3 (dropping quadratic
#terms).
interquad <- function(x, z){
  #define the base model to which the reduced models will be compared
  base <- x
  #defining the z vector without "mes", bc "mes" does not have interactions
  #nor quadratic terms
  z <- z[z != "Mes"]
  #creating a vector with the names of all the interactions
  int <- attr(base$terms, "term.labels")[grepl(":", attr(base$terms, "term.labels"))]
  #Create variable w_com that serves as the condition for the while-loop
  #After the first run, it is equal to the "best" model
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #less complex model is better. Second condition makes sure that there are 
  #actually interaction terms.
  while((rownames(w_com)[1] != "base") & (length(int) > 0)){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(int) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in int){
      #index
      j <-match(i,int)
      #We run the model without the interaction terms of i and store the 
      #fitted model in st
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(i, attr(base$terms, "term.labels")),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(int, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs, The first row contains the preferred model
    #(the one with the lowest WAIC)
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The interaction terms of variable ", rownames(w_com)[1], " will be dropped"))
    
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[1])]]
    #deleting the dropped variable in the names vector (Es sind ja immer zwei,
    #weil ja zwei Variablen mit einem Interaction term assoziert sind)
    int <- int[int != rownames(w_com)[1]] 
  }
  #Now doing more or less the same with the quadratic terms
  #resetting w_com variable
  w_com <- matrix("start")
  rownames(w_com) <- "start"
  #initializing a while-loop that keeps dropping interactions as long as the 
  #one less complex model is better
  while(rownames(w_com)[1] != "base"){
    #Create an empty list to store the model-objects of the different models
    st <- vector("list", length(z) + 1)
    #an empty list for all the WAICs
    waics <- st
    #add the baseline WAIC to which the more parsimonious models will 
    #be compared (It's always the "Most" complex model at the specific point of the 
    #model-selection procedure)
    waics[[length(waics)]] <- waic(base)
    #Calculate the different reduced models
    for (i in z){
      #index
      j <-match(i,z)
      #We run the model without the terms associated with z (linear, quadratic
      #and interaction terms and store the fitted model in st)
      st[[j]] <- update(base, formula=drop.terms(base$terms, grep(paste0("I(", i), attr(base$terms, "term.labels"), fixed = T),
                                                 keep.response=TRUE))
      print(attr(st[[j]]$terms, "term.labels"))
      #Calculating the WAIC of the model and storing it in waics
      waics[[j]] <- waic(st[[j]])
    }
    #a vector with all the model names (the ones without one covariate + the base)
    n <- c(z, "base")
    #naming the two storage lists accordingly
    names(st) <- n
    names(waics) <- n
    #Compare the different WAICs, The first row contains the preferred model
    #(the one with the lowest WAIC)
    w_com <- loo_compare(waics)
    print(w_com)
    #Printing a message to the screen indicating that variable xxx will be thrown out
    print(paste0 ("The quadratic term of variable ", rownames(w_com)[1], " will be dropped"))
    
    #storing the model without the removed variable as the new base model
    base <- st[[which(names(st) == rownames(w_com)[1])]]
    
    #deleting the variable that we throw out according to WAIC in the variable
    #names vector z
    z <- z[z != rownames(w_com)[1]]
  }
  #Return the final model
  base#for some reason this does not work
}