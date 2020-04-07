
# generic methods for "GPA" class
setGeneric(
  "get_fit",
  function(x) standardGeneric("get_fit")
)

setGeneric(
  "get_setting",
  function(x) standardGeneric("get_setting")
)

setGeneric(
  "get_gwasPval",
  function(x) standardGeneric("get_gwasPval")
)

setGeneric(
  "get_annMat",
  function(x) standardGeneric("get_annMat")
)

setGeneric( "fdr",
    function( object, ... )
    standardGeneric("fdr")
)

setGeneric( "assoc",
    function( object, ... )
    standardGeneric("assoc")
)

setGeneric( "cov",
    function( object, ... )
    standardGeneric("cov")
)

setGeneric( "se",
    function( object, ... )
    standardGeneric("se")
)

setGeneric( "estimates",
    function( object, ... )
    standardGeneric("estimates")
)
