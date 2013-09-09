# Note: Curated models can be saved as an R-list via get_all_models(),
# but this was not possible with non-curated models. (the '_noncur'
# functions extract only the required information required from each model).

# Function get_all_models():
# Reads a file containing SBML models from the curated section of the
# BioModels database and returns these models in an R list
get_all_models = function() {
  l = list.files("./BioModels/curated/")
  lmod = list()
    
  for(i in 1:length(l)) {
    filename = paste0("./BioModels/curated/", l[[i]])
    # no need to set defaults for analysis, so use libSBML's readSBML 
    doc = libSBML::readSBML(filename)
    lmod[[i]] = SBMLDocument_getModel(doc)
  }
  
  return(lmod)
}



# ================================================================
# Model history
# ================================================================

# Function get_history_cur(lmod):
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a character vector containing creation dates (where available)
# or message string indicating that history / creation date not set
get_history_cur = function(lmod) {
  
  num = length(lmod)
  v = vector(mode = "character")
  
  for(i in 1:num) {
    date = get_creation_data(lmod[[i]])
    v = append(v, date)
  }
  
  return(v)
}



# Function get_history_noncur():
# Reads a file containing SBML models from the noncurated section of the
# BioModels database and checks model history
# Returns a character vector containing creation dates (where available)
# or message string indicating that history / creation date not set
get_history_noncur = function() {
  
  l = list.files("./BioModels/noncurated/")
  v = vector(mode = "character")
  
  for(i in 1:length(l)) {
    filename = paste0("./BioModels/noncurated/", l[[i]]) 
    doc = libSBML::readSBML(filename)
    mod = SBMLDocument_getModel(doc)    
    date = get_creation_data(mod)
    v = append(v, date)    
    delete(doc)
  }
  
  return(v)
}



# Function get_creation_data(mod):
# Determines whether or not a creation date is set for the model
# Argument 'mod' is a model from the BioModels database
# Returns model creation date as a string if set,
# otherwise returns a message string
get_creation_data = function(mod) {
  
  if(SBase_isSetModelHistory(mod)) {
    mh = SBase_getModelHistory(mod)
  } else return ("History not set")
  
  if(ModelHistory_isSetCreatedDate(mh)) {
    cd = ModelHistory_getCreatedDate(mh)
    date = Date_getDateAsString(cd)
    # trim date string to YYYY-MM-DD format
    return(strsplit(date, "T")[[1]][1])    
  } else return ("Creation Date not set")
}



# Function track_num_models(v):
# Tracks cumulative number of models using the model creation date
# Argument v is a character vector containing model creation dates
# Returns a list of dates with the corresponding cumulative sum
# N.B. Not all models have creation dates
track_num_models = function(v) {
  
  # sort contents of v as dates
  v = sort(as.Date(v))
  total = vector(mode = "numeric")
  
  for(i in 1:length(v)) {
    total[i] = sum(v <= v[i])
  }
  
  return(list(v,total))
}



# ==================================================================
# Level and Version data
# ==================================================================

# Function get_lvd_cur(lmod):
# Get level, version and creation date of curated models
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a data frame containing level, version and creation date
# (where available) for each model in the list 'lmod'
get_lvd_cur = function(lmod) {
  
  num = length(lmod)
  df = data.frame(Level= integer(0),
                  Version= integer(0),
                  Date = character(0))  
  
  for(i in 1:num) {
    df = rbind(df, get_level_version_date(lmod[[i]]))    
  }
  return(df)
}



# Function get_lvd_noncur(lmod):
# Get level, version and creation date of noncurated models
# Returns a data frame containing level, version and creation date
# (where available)
get_lvd_noncur = function() {
  
  l = list.files("./BioModels/noncurated/")
  df = data.frame(Level= integer(0),
                  Version= integer(0),
                  Date = character(0))  
  
  for(i in 1:length(l)) {
    filename = paste0("./BioModels/noncurated/", l[[i]]) 
    doc = libSBML::readSBML(filename)
    mod = SBMLDocument_getModel(doc)    
    df = rbind(df, get_level_version_date(mod))
    delete(doc)
  }
  
  return(df)
}


# Function get_level_version_date(mod):
# Gets level, version and creation date of a model
# Argument 'mod' is a models from the BioModels database
# Returns a data frame containing level, version and creation date
get_level_version_date = function(mod) {
  
  level = SBase_getLevel(mod)
  vers = SBase_getVersion(mod)    
  date = get_creation_data(mod)
  return(data.frame(level,vers,date))
}



# ==================================================================
# Distribution of Components
# ==================================================================

# Function get_component_distribution_cur(lmod):
# Gets the number of occurrences of each component in all curated models
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a matrix with rows corresponding to each model in the list,
# and columns corresponding to each component
get_component_distribution_cur = function(lmod) {
  
  num = length(lmod)
  x = vector(mode = "numeric")
  
  for(i in 1:num) {
    # use rbind to build matrix
    x = rbind(x, get_num_components(lmod[[i]]))
  }
  
  return(x)
}



# Function get_component_distribution_noncur():
# Gets the number of occurrences of each component in all noncurated models
# Returns a matrix with rows corresponding to each model in the list,
# and columns corresponding to each component
get_component_distribution_noncur = function() {
  
  l = list.files("./BioModels/noncurated/")
  x = vector(mode = "numeric")
  
  for(i in 1:length(l)) {    
    filename = paste0("./BioModels/noncurated/", l[[i]]) 
    doc = libSBML::readSBML(filename)
    mod = SBMLDocument_getModel(doc)    
    # use rbind to build matrix
    x = rbind(x, get_num_components(mod))   
    delete(doc)
  }
  
  return(x)  
}



# Function get_num_components(mod):
# Gets the number of occurrences of each component in a model
# Argument 'mod' is a model from the BioModels database
# Returns a vector of named elements
get_num_components = function(mod) {
  
  components = c("FunctionDefinitions", "UnitDefinitions", "CompartmentTypes",
                 "SpeciesTypes", "Compartments", "Species", "Parameters",
                 "InitialAssignments", "Rules", "Constraints",
                 "Reactions", "Events")  
  
  v = vector(mode = "numeric")
  
  for(comp in components) {
    v[comp] = do.call(paste0("Model_getNum", comp), list(mod))
  }
  
  return(v)
}



# ==================================================================
# Distribution of SubComponents
# (components within the Unit, Reaction and Event components)
# ==================================================================

# Function get_subComponent_distribution_cur(lmod):
# Gets the number of occurrences of subcomponents (ie components
# within components) in all curated models
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a matrix with rows corresponding to each model in the list,
# and columns corresponding to each subcomponent
get_subComponent_distribution_cur = function(lmod) {
  
  num = length(lmod)
  x = vector(mode = "numeric")
  
  for(i in 1:num) {
    # use rbind to build matrix
    x = rbind(x, get_model_subComponents(lmod[[i]]))
  }
  
  return(x)
}



# Function get_model_subComponents(mod):
# Gets the number of occurrences of each subcomponent in a model
# Argument 'mod' is a model from the BioModels database
# Returns a vector of named elements
get_model_subComponents = function(mod) {
  
  return(c(get_reaction_components(mod),
           get_event_components(mod),
           get_unitDefinition_components(mod)))
}



# Function get_reaction_components(mod)
# Gets number of Reaction sub-components (Modifiers, Products,
# Reactants, KineticLaws, LocalParameters) in a model
get_reaction_components = function(mod) {
  
  components = c("Modifiers", "Products", "Reactants")  
  v = vector(mode = "numeric")
  numReactions = Model_getNumReactions(mod)
  # variables to store cumulative sums
  numModifiers = 0
  numProducts = 0
  numReactants = 0
  numKineticLaws = 0
  numLocalParameters = 0
  
  if(numReactions > 0) {
    for(j in 1:numReactions) {
      react = mod[[11]][[j]]
      for(comp in components) {
        v[comp] = do.call(paste0("Reaction_getNum", comp), list(react))
      }
      # get cumulative sums for each subcomponent
      numModifiers = numModifiers + v[["Modifiers"]]
      numProducts = numProducts + v[["Products"]]
      numReactants = numReactants + v[["Reactants"]]
      # now check for KineticLaw
      if(Reaction_isSetKineticLaw(react)) {
        kl = Reaction_getKineticLaw(react)
        numKineticLaws = numKineticLaws + 1
        numLocalParameters = numLocalParameters + KineticLaw_getNumLocalParameters(kl)            
      }      
    }       
  }
  return(c("Modifier" = numModifiers,
           "Product" = numProducts,
           "Reactant" = numReactants,
           "KineticLaw" = numKineticLaws,
           "LocalParameter" = numLocalParameters))
}



# Function get_event_components(mod)
# Gets number of Event sub-components (Trigger, Delay, Prioirity,
# EventAssignments) in a model
get_event_components = function(mod) {
  
  numEvents = Model_getNumEvents(mod)
  # variables to store cumulative sums
  numTriggers = 0
  numDelays = 0
  numPriority = 0
  numEventAssignments = 0
  
  if(numEvents > 0) {
    for(j in 1:numEvents) {
      event = mod[[12]][[j]]
      if(Event_isSetTrigger(event)) {
        numTriggers = numTriggers + 1
      }
      if(Event_isSetDelay(event)) {
        numDelays = numDelays + 1
      }
      if(Event_isSetPriority(event)) {
        numPriority = numPriority + 1
      }
      numEventAssignments = numEventAssignments + Event_getNumEventAssignments(event)
    }    
  }
  
  return(c("Trigger" = numTriggers,
           "Delay" = numDelays,
           "Priority" = numPriority,
           "EventAssignment" = numEventAssignments))
}



# Function get_unitDefinition_components(mod)
# Gets number of UnitDefinition sub-components (Units) in a model
get_unitDefinition_components = function(mod) {
  
  numUnitDefs = Model_getNumUnitDefinitions(mod)
  numUnits = 0
  
  if(numUnitDefs > 0) {
    for(j in 1:numUnitDefs) {
      numUnits = numUnits + UnitDefinition_getNumUnits(mod[[2]][[j]])      
    }    
  }
  
  return(c("Unit" = numUnits))
}



# ==================================================================
# Connections (Species)
# ==================================================================

# Function get_species_resources_cur(lmod):
# Gets the URIs for all species in all curated models
# N.B. only URIs for which the biological qualifier is 
# BQB_IS or BQB_IS_VERSION_OF are included
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a vector of species URIs
# Vector elements are named with the id of the model
# in which the URIs exist
get_species_resources_cur = function(lmod) {
  
  num = length(lmod)
  v = vector(mode = "character")
  
  for(i in 1:num) {
    v = append(v, get_model_species_resources(lmod[[i]]))
  }
  
  return(v)  
}


# returns a vector containing species URIs for a model (duplicates 
# within a model are disregarded)
# The model id is used to name the elements of the vector 
get_model_species_resources = function(mod) {
  
  v = vector(mode = "character")
  
  numSpecies = Model_getNumSpecies(mod)
  if(numSpecies > 0) {
    for(i in 1:numSpecies) {
      v = append(v, get_species_resources(mod[[6]][[i]]))    
    }    
  }  
  # disregard duplicates within a model
  v = unique(v)
  
  # name elements and return
  if(length(v)>0) {
    names(v)[1:length(v)] = Model_getId(mod)
  }
  return(v)  
}

# returns a vector containing all URIs for a species for which
# the biological qualifier is BQB_IS or BQB_IS_VERSION_OF
get_species_resources = function(sp) {
  
  v = vector(mode = "character")
  
  nCVT = SBase_getNumCVTerms(sp)
  
  if(nCVT > 0) {    
    
    for(i in 1:nCVT) {
      
      cvt = SBase_getCVTerm(sp,i-1)
      
      if(CVTerm_hasRequiredAttributes(cvt)){
        
        bqType = CVTerm_getBiologicalQualifierType(cvt)
        
        # only add to vector if qualifier is BQB_IS or BQB_IS_VERSION_OF
        if(bqType %in% c("BQB_IS", "BQB_IS_VERSION_OF")) {
          nRes = CVTerm_getNumResources(cvt)
          for(j in 1:nRes) {
            v = append(v, CVTerm_getResourceURI(cvt,j-1))
          }          
        }             
      }
    }    
  }  
  return(v)
}


# ==================================================================
# Connections (Reactions)
# ==================================================================

# Function get_reaction_resources_cur(lmod):
# Gets the URIs for all reactions in all curated models
# N.B. only URIs for which the biologival qualifier is 
# BQB_IS or BQB_IS_VERSION_OF are included
# Argument 'lmod' is list of curated models from the BioModels database
# Returns a vector of species URIs
# Vector elements are named with the id of the model
# in which the URIs exist
get_reaction_resources_cur = function(lmod) {
  
  num = length(lmod)
  v = vector(mode = "character")
  
  for(i in 1:num) {    
    v = append(v, get_model_reaction_resources(lmod[[i]]))
  }
  
  return(v)
}



# returns a vector containing resource URIs for a model (duplicates 
# within a model are disregarded)
# The model id is used to name the elements of the vector 
get_model_reaction_resources = function(mod) {
  
  v = vector(mode = "character")
  
  numReactions = Model_getNumReactions(mod)
  if(numReactions > 0) {
    for(i in 1:numReactions) {
      v = append(v, get_reaction_resources(mod[[11]][[i]]))    
    }    
  }  
  # disregard duplicates within a model
  v = unique(v)
  
  # name elements and return
  if(length(v)>0) {
    names(v)[1:length(v)] = Model_getId(mod)
  }
  return(v)  
}



# returns a vector containing all URIs for a reaction for which
# the biological qualifier is BQB_IS or BQB_IS_VERSION_OF
get_reaction_resources = function(r) {
  
  v = vector(mode = "character")
  
  nCVT = SBase_getNumCVTerms(r)
  
  if(nCVT > 0) {    
    
    for(i in 1:nCVT) {
      
      cvt = SBase_getCVTerm(r,i-1)
      
      if(CVTerm_hasRequiredAttributes(cvt)){
        
        bqType = CVTerm_getBiologicalQualifierType(cvt)
        
        # only add to vector if qualifier is BQB_IS or BQB_IS_VERSION_OF
        if(bqType %in% c("BQB_IS", "BQB_IS_VERSION_OF")) {
          nRes = CVTerm_getNumResources(cvt)
          for(j in 1:nRes) {
            v = append(v, CVTerm_getResourceURI(cvt,j-1))
          }          
        }             
      }
    }    
  }  
  return(v)
}
  

# ==================================================================
# Connections (URIs occurring in more than one model)
# ==================================================================

# Function get_component_resources_cur(lmod, compStr):
# Gets the URIs for particular components in all curated models
# N.B. only URIs for which the biological qualifier is 
# BQB_IS or BQB_IS_VERSION_OF are included
# Argument 'lmod' is list of curated models from the BioModels database
# Argument 'compStr' is an SBML model component name (plural)
# Returns a vector of component URIs
# Vector elements are named with the id of the model
# in which the URIs exist
get_component_resources_cur = function(lmod, compStr) {
  
  componentNames = c("FunctionDefinitions", "UnitDefinitions",
                     "CompartmentTypes", "SpeciesTypes", "Compartments",
                     "Species", "Parameters", "InitialAssignments",
                     "Rules", "Constraints", "Reactions", "Events") 
  
  if(!(compStr %in% componentNames)) {
    messageStr = paste0('"', compStr, '" is not a valid component string.')
    stop(messageStr)    
  }
    
  num = length(lmod)
  v = vector(mode = "character")
  
  for(i in 1:num) {
    
    cat("\nAt model:", i)
    
    v = append(v, get_model_component_resources(lmod[[i]], compStr))
  }
  
  return(v)
}


# returns a vector containing species URIs for a model (duplicates 
# within a model are disregarded)
# The model id is used to name the elements of the vector 
get_model_component_resources = function(mod, compStr) {
  
  v = vector(mode = "character")
    
  numComponent = do.call(paste0("Model_getNum", compStr), list(mod))
  if(numComponent > 0) {
    for(i in 1:numComponent) {
      v = append(v, get_component_resources(mod[[paste0("ListOf",compStr)]][[i]]))    
    }    
  }  
  # disregard duplicates within a model
  v = unique(v)
  
  # name elements and return
  if(length(v)>0) {
    names(v)[1:length(v)] = Model_getId(mod)
  }
  return(v)  
}

# returns a vector containing all URIs for a species for which
# the biological qualifier is BQB_IS or BQB_IS_VERSION_OF
get_component_resources = function(obj) {
  
  v = vector(mode = "character")
  
  nCVT = SBase_getNumCVTerms(obj)
  
  if(nCVT > 0) {    
    
    for(i in 1:nCVT) {
      
      cvt = SBase_getCVTerm(obj,i-1)
      
      if(CVTerm_hasRequiredAttributes(cvt)){
        
        bqType = CVTerm_getBiologicalQualifierType(cvt)
        
        # only add to vector if qualifier is BQB_IS or BQB_IS_VERSION_OF
        if(bqType %in% c("BQB_IS", "BQB_IS_VERSION_OF")) {
          nRes = CVTerm_getNumResources(cvt)
          for(j in 1:nRes) {
            v = append(v, CVTerm_getResourceURI(cvt,j-1))
          }          
        }             
      }
    }    
  }  
  return(v)
}
