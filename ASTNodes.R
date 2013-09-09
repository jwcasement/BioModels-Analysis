# These functions walk an AST and determine each node type

# use a global to store AST types
types = vector(mode = "character")

# define constants and basic functions to be disregarded
disregard = c("AST_PLUS", "AST_MINUS", "AST_TIMES", "AST_DIVIDE",
              "AST_REAL", "AST_INTEGER", "AST_RATIONAL", "AST_NAME")

# get AST types for all nodes with 'astNode' as a root
get_all_types = function(astNode) {
  if(length(types) > 0) types <<- vector(mode = "character")
  walk_tree(astNode)
  # remove anything in the 'disregard' vector, and any duplicates
  return(unique(types[!(types %in% disregard)]))
}

# a recursive function to walk the AST and get node types
# node types are added to global vector 'types' 
walk_tree = function(astNode) {
  
  n = ASTNode_getNumChildren(astNode)
  types[length(types)+1] <<- ASTNode_getType(astNode)
    
  if(n > 0) {    
    for(i in 1:n) {      
      child = ASTNode_getChild(astNode, i-1)
      walk_tree(child)      
    }
  }
}
