# These functions determine AST node types in a list of models

library(libsbmlwrapper)
source("./BiomodelsAnalysis/ASTNodes.R")

# get all AST nodes types in a list of models
get_AST_from_list = function(lmod, removeModelDuplicates) {
    
  num = length(lmod)
  allTypes = vector(mode = "character")
  
  for(i in 1:num) {
    allTypes = append(allTypes, get_all_AST(lmod[[i]], removeModelDuplicates))    
  }
    
  return(allTypes)
}



# get all AST node types in a model
# if removeModelDuplicates is TRUE, replication of an AST node type within
# a model is disregarded
get_all_AST = function(mod, removeModelDuplicates) {
  
  v = vector(mode = "character")
  
  # loop through all components which may contain Math elements
  for(compStr in c("FunctionDefinition", "InitialAssignment",
                   "Rule", "Constraint", "Reaction", "Event")) {
    # check that component exists in model
    size = do.call(paste0("Model_getNum", compStr, "s"), list(mod))
    if(size > 0) {      
      for(i in 1:size) {
        # get the AST types for the component object and update vector
        v = append(v, get_AST_for_component(compStr, mod, i))        
      }
    }
  }
  ifelse(removeModelDuplicates == TRUE, return(unique(v)), return(v))
}


# Function returns a vector of all AST node types (duplicates removed)
# from the Math element within a component 
get_AST_for_component = function(compStr, model, i) {
  
  v = vector(mode = "character")
  
  if(compStr == "Reaction") {
    # for Reaction component, Math element exists within KineticLaw component
    r = model[[11]][[i]]
    # get KineticLaw object
    if(Reaction_isSetKineticLaw(r)) {
      kl = r[["KineticLaw"]]
      # check that Math element is set for KineticLaw object
      if(KineticLaw_isSetMath(kl)) {
        # if Math is set, update vector of AST nodes
        AST = get_all_types(kl[["Math"]])
        if(length(AST)>0) {
          names(AST)[1:length(AST)] = "KineticLaw"
          v = append(v, AST)
        }
      }
    }
  } else if(compStr == "Event") {
    # for Event component, Math element exists within Trigger, Delay,
    # Priority or EventAssignment subcomponents
    ev = model[[12]][[i]]
    # check for Trigger objects
    if(Event_isSetTrigger(ev)) {
      tr = ev[["Trigger"]]
      # check that Math element is set for Trigger object
      if(Trigger_isSetMath(tr)) {
        # if Math is set, update vector of AST nodes
        AST = get_all_types(tr[["Math"]])
        if(length(AST)>0) {
          names(AST)[1:length(AST)] = "Trigger"
          v = append(v, AST)
        }                
      }
    }
    # check for Delay objects
    if(Event_isSetDelay(ev)) {
      delay = ev[["Delay"]]
      # check that Math element is set for Delay object
      if(Delay_isSetMath(delay)) {
        # if Math is set, update vector of AST nodes
        AST = get_all_types(delay[["Math"]])
        if(length(AST)>0) {
          names(AST)[1:length(AST)] = "Delay"
          v = append(v, AST)
        }                
      }      
    }
    # check for Priority objects
    if(Event_isSetPriority(ev)) {
      pr = ev[["Priority"]]
      # check that Math element is set for Priority object
      if(Priority_isSetMath(pr)) {
        # if Math is set, update vector of AST nodes
        AST = get_all_types(pr[["Math"]])
        if(length(AST)>0) {
          names(AST)[1:length(AST)] = "Priority"
          v = append(v, AST)
        }        
      }      
    }
    # check for EventAssignments
    size = Event_getNumEventAssignments(ev)
    if(size > 0) {
      for(i in 1:size) {
        # get EventAssignment object
        ea = ev[["EventAssignment"]][[i]]
        if(EventAssignment_isSetMath(ea)) {
          # if Math is set, update vector of AST nodes
          AST = get_all_types(ea[["Math"]])
          if(length(AST)>0) {
            names(AST)[1:length(AST)] = "EventAssignment"
            v = append(v, AST)
          }        
        }
      }
    }    
  } else {
    # for other components, Math element exists directly within component itself 
    comp = model[[paste0("ListOf", compStr, "s")]][[i]]
    # check that Math element is set for component object
    if(do.call(paste0(compStr, "_isSetMath"), list(comp))) {
      # if Math is set, update vector of AST nodes
      AST = get_all_types(comp[["Math"]])
      if(length(AST) > 0) {
        names(AST)[1:length(AST)] = compStr
        v = append(v, AST)
      }      
    }
  }
  # remove duplicates (keeping names) and return
  return(v[!duplicated(v)])
}
