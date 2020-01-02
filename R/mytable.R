#' @title mytable
#' @name mytable
#' @description build a contingency table of the counts at each combination of factor levels with two objects
#' @param x one object
#' @param y one object
#' @return a table \code{n}
#' @examples
#' \dontrun{
#' x<-c(4,5,6,7,8)
#' y<-c(1,2,3,3,4)
#' mytable(x,y)
#' }
#' @export
mytable <-function(x,y){
  x=as.factor(x)
  y=as.factor(y)
  xi=as.integer(x)
  xm=as.integer(max(xi))
  yi=(as.integer(y)-1L)*xm
  ym=as.integer(max(yi)+xm)
  matrix(.Internal(tabulate(xi+yi,ym)),xm)
}
