\name{wrapper functions}
\alias{get.fn}
\alias{get.gr}
\alias{get.fngr}
\alias{get.hessian}
\title{Compute function, gradient and Hessian from objects of class sparseHessianObj.}
\description{Wrapper functions to call methods of an object of class sparseHessianObj.}
\usage{
get.fn(x, obj)
get.gr(x, obj)
get.fngr(x, obj)
get.hessian(x, obj)
}
\arguments{
	\item{x}{Value at which fn, gr or the hessian is evaluated.}
	\item{obj}{Object of class sparseHessianObj that contains all of the function, gradient and Hessian information.}
}
