collapse_by_taxonomy.default <- function(Tab,Tax,Group = NULL, level=4,
                                         FUN=sum,sepchar=";"){
  # Match taxonomy and table rows
  #row.names(Tax) <- as.character(Tax$ID)
  
  #if(class(Tab) != "matrix")
    #stop("ERROR: Tab must be a matrix object",call.=TRUE)
  if(nrow(Tab) != nrow(Tax))
    stop("ERROR: Number of rows in Tab and Tax do not match",call.=TRUE)
  if(any(row.names(Tab) != row.names(Tax)))
    stop("ERROR: Row names for Tab and Tax do not match",call.=TRUE)
  
  if(is.null(Group) && is.numeric(level)){
    tax <- get_tax_level(Tax,level=level,sepchar=sepchar)
  }else if(!is.null(Group)){
    if(length(Group) == 1){
      tax <- Tax[ ,Group ]
    }else{
      tax <- factor(Group)
    }
  }else{
    stop("ERROR: Either a non-null group must be provided, or a numeric level for taxonomyu collapsint",
         call. = TRUE)
  }
  
  
  Tab.collapsed <- collapse_matrix(x=Tab,groups=tax,dim=1,FUN=FUN)
  return(Tab.collapsed)
}


#' @rdname collapse_by_taxonomy
#' @method collapse_by_taxonomy Dataset
#' @export
collapse_by_taxonomy.Dataset <- function(Dat, Group = NULL, level = 4,
                                         FUN = sum, sepchar = ";"){
  res <- collapse_by_taxonomy.default(Tab = Dat$Tab,
                                      Tax = Dat$Tax,
                                      Group = Group,
                                      level = level,
                                      FUN = FUN,
                                      sepchar = sepchar)
  
  if(length(Dat$Map) > 0){
    res <- create_dataset(Tab = res, Map = Dat$Map)
  }
  
  return(res)
}
