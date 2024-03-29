% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tan.dependence.Diffuse_S.mult.R
\name{Tan.dependence.Diffuse_S.mult}
\alias{Tan.dependence.Diffuse_S.mult}
\title{计算林下多个位点的加权进界邻体散射荫蔽度}
\usage{
Tan.dependence.Diffuse_S.mult(a, b, tan, MI, wei, tan2)
}
\arguments{
\item{a}{林下研究位点的横坐标、纵坐标与高度}

\item{b}{样地内林木的横坐标与树高}

\item{tan}{局域临界值}

\item{MI}{校正系数}

\item{wei}{线性校正系数}

\item{tan2}{背景临界值，可以省略}
}
\value{
林下多个位点的加权进界邻体散射荫蔽度
}
\description{
计算林下多个位点的加权进界邻体散射荫蔽度
}
\examples{
a=data.frame(x=c(25,35,40),y=c(25,21,10),size=c(0,1,2))
b=Understory.Light::b
Understory.light::Tan.dependence.Diffuse_S.mult(a=a,b=b,tan=2.8,MI=10.660,wei=0.5)
}
