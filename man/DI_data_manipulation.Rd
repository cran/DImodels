\name{DI_data_manipulation}
\alias{DI_data_manipulation}
\alias{DI_data_ADD}
\alias{DI_data_E_AV}
\alias{DI_data_FG}
\alias{DI_data_prepare}
\alias{DI_data_fullpairwise}
\encoding{UTF-8}
\title{
Data manipulation functions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The following five functions compute additional variables for the various types of interactions among pairs of species proportions. These variables are defined following Kirwan et al 2007 and 2009, and Connolly et al 2013. 
 
\code{DI_data_E_AV}: creates the average pairwise interaction variable (AV) and a scaled version of it (E). The variable AV is the sum of products of the  proportions of each pair of species in the mixture. The variable E is a scaled version of AV that ranges between 0 (for a monoculture community) to 1 for the equi-proportional mixture of all species in the pool.  
 
\code{DI_data_FG}: creates the functional group interaction variables. There is a variable for (within) each functional group and one for (between) each pair of functional groups, i.e., if there are two functional groups, there will be three functional group interaction variables, while if there are three functional groups, there will be six functional group interaction variables.  
 
\code{DI_data_ADD}: creates the additive species interaction variables, one for each species. 

\code{DI_data_prepare}: computes all types of interactions in the preceding three functions. Use this function to implement all three previous functions in one step (avoiding the need to use the individual ones).

\code{DI_data_fullpairwise}: computes all individual pairwise interactions. There will be \eqn{s*(s-1)/2} new variables created, where \emph{s} is the number of species in the pool.


By default, the interaction variables described above are created with \code{theta = 1}, but a different value of theta can also be specified (Connolly et al 2013). 
}
\usage{
DI_data_E_AV(prop, data, theta = 1)
DI_data_FG(prop, FG, data, theta = 1)
DI_data_ADD(prop, data, theta = 1)
DI_data_prepare(y, block, density, prop, treat, FG = NULL, data, theta = 1)
DI_data_fullpairwise(prop, data, theta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prop}{ 
%%     ~~Describe \code{P} here~~
A vector of column names identifying the species proportions in the dataset. For example, if the species proportions columns are labelled p1 to p4, then \code{prop = c("p1","p2","p3","p4")}. The column numbers in which the proportions are stored can also be referred to, for example,  \code{prop = 4:7} for the \code{Switzerland} data.
}
  \item{FG}{
%%     ~~Describe \code{FG} here~~
If species are classified by \emph{g} functional groups, this argument takes a text list (of length \emph{s}) of the functional group to which each species belongs. For example, for four grassland species with two grasses and two legumes, it could be \code{FG = c("G","G","L","L")}, where G stands for grass and L stands for legume. This argument is required in \code{DI_data_FG}. This argument is not required in \code{DI_data_prepare}, but if omitted, the functional group interactions will not be computed by the function. 
}
  \item{data}{
%%     ~~Describe \code{data} here~~
Specify the dataset, for example, \code{data = Switzerland}. The dataset name should not appear in quotes. 
}
\item{theta}{
%%     ~~Describe \code{data} here~~
Interaction variables will be computed with the theta power, equal to the value specified, on all \eqn{pi*pj} components of each interaction variable, with default value one. For example, with three species \code{AV = p1*p2 + p1*p3 + p2*p3} and if computed with \code{theta} = 0.5, this becomes \code{(p1*p2)^0.5 + (p1*p3)^0.5 + (p2*p3)^0.5}.
}
\item{y}{
%%     ~~Describe \code{y} here~~
The column name of the response vector, for example, \code{y = "yield"}. The name of the response variable must be contained in quotes. This argument is not required (but is used internally in the \code{DI} and \code{autoDI} functions). If this argument is omitted, there are no implications for the calculation of interaction variables here.
}
\item{block}{
%%     ~~Describe \code{block} here~~
The column name of the block variable. This argument is not required (but is used internally in the \code{DI} and \code{autoDI} functions). If there is no blocking variable, omit this argument; in this case, a new variable \code{block_zero}, which is a column of zeros, will be computed. This column of zeros is used internally in the \code{DI} and \code{autoDI} functions, but has no implications the calculation of interaction variables here. 
}
\item{density}{
%%     ~~Describe \code{block} here~~
The column name of the density variable. This argument is not required (but is used internally in the \code{DI} and \code{autoDI} functions). If there is no density variable, omit this argument; in this case, a new variable \code{density_zero}, which is a column of zeros, will be computed. This column of zeros is used internally in the \code{DI} and \code{autoDI} functions, but has no implications for the calculation of interaction variables here. 
}
\item{treat}{
%%     ~~Describe \code{treat} here~~
The column name of the treatment variable. This argument is not required (but is used internally in the \code{DI} and \code{autoDI} functions). If there is no treatment variable, omit this argument; in this case, a new variable \code{treat_zero}, which is a column of zeros, will be computed. This column of zeros is used internally in the \code{DI} and \code{autoDI} functions, but has no implications for the calculation of interaction variables here. 
}

}
\details{ 
%%  ~~ If necessary, more details than the description above ~~

 
\strong{What are Diversity-Interactions models?} 
 

Diversity-Interactions (DI) models (Kirwan et al 2009) are a set of tools for analysing and interpreting data from experiments that explore the effects of species diversity on community-level responses. We recommend that users of the \code{DImodels} package read the short introduction to DI models (available at: \code{\link{DImodels}}). Further information on DI models is available in Kirwan et al 2009 and Connolly et al 2013.

 
\strong{Checks on data prior to using the data manipulation functions.} 
 

Before applying the data manipulation functions to your dataset, check that the species proportions in each row sum to one. See the 'Examples' section for code to do this. An error message will be generated if the proportions don't sum to one. 

 
\strong{When are the data manipulation functions needed?} 
 

It is not required to use the data manipulation functions if using the \code{autoDI} function, or the \code{DImodel} option in the \code{DI} function, as they will automatically create the species interaction variables needed. If using species interaction variables in the \code{extra_formula} or \code{custom_formula} options in \code{DI}, then it is required to have the variables already in the dataset and these functions can do that.  

 
\strong{Short worked example to illustrate how the data manipulation functions work} 
 

The code to implement this example is provided in the 'Examples' section. 

Assume four species with initial proportions in two communities: (0.1, 0.2, 0.3, 0.4) and (0.25, 0.25, 0.25, 0.25), with \code{FG = c("G","G","L","L")}. 
 
For community 1: (0.1,0.2,0.3,0.4), assuming theta = 1, the data preparation functions will compute the following additional variables (details in Kirwan et al 2007 and 2009): 

AV = 0.1*0.2 + 0.1*0.3 + 0.1*0.4 + 0.2*0.3 + 0.2*0.4 + 0.3*0.4 = 0.35

E = \emph{(2s/(s-1))}*AV = 0.9333

p1_add = 0.1 * (1 - 0.1) = 0.09

p2_add = 0.2 * (1 - 0.2) = 0.16

p3_add = 0.3 * (1 - 0.3) = 0.21

p4_add = 0.4 * (1 - 0.4) = 0.24

bfg_G_L = 0.1*0.3 + 0.1*0.4 + 0.2*0.3 + 0.2*0.4 = 0.21

wfg_G = 0.1*0.2 = 0.02

wfg_L = 0.3*0.4 = 0.12

For community 1: (0.1,0.2,0.3,0.4), assuming theta = 0.5, the data preparation functions will compute the follow additional variables (details in Connolly et al 2013): 

AV = (0.1*0.2)^0.5 + (0.1*0.3)^0.5 + (0.1*0.4)^0.5 + (0.2*0.3)^0.5 + (0.2*0.4)^0.5 + (0.3*0.4)^0.5 = 1.3888
 
E =\emph{(2s/(s-1))}*AV = 3.7035

p1_add = 0.1^0.5 * (0.2^0.5 + 0.3^0.5 + 0.4^0.5) = 0.5146

p2_add = 0.2^0.5 * (0.1^0.5 + 0.3^0.5 + 0.4^0.5) = 0.6692

p3_add = 0.3^0.5 * (0.1^0.5 + 0.2^0.5 + 0.4^0.5) = 0.7646

p4_add = 0.4^0.5 * (0.1^0.5 + 0.2^0.5 + 0.3^0.5) = 0.8293

bfg_G_L = (0.1*0.3)^0.5 + (0.1*0.4)^0.5 + (0.2*0.3)^0.5 + (0.2*0.4)^0.5 = 0.9010

wfg_G = (0.1*0.2)^0.5 = 0.1414

wfg_L = (0.3*0.4)^0.5 = 0.3464  

When using the data manipulation functions to create interactions for theta values for a value different from 1, it is recommended to rename the new interaction variables to include \code{_theta}. 

The data manipulation values for community 2 can be seen when the 'Examples' section code is run.
}
\value{
The \code{DI_data_prepare} function returns a named list with the following components:
\item{newdata}{a \code{data.frame} containing all manipulated variables}
\item{y}{the response variable name}
\item{block}{the block variable name}
\item{density}{the density variable name}
\item{prop}{the species proportions variable names}
\item{treat}{the treatment variable name}
\item{FG}{the variables used in the FG model}
\item{P_int_flag}{a logical value used internally}
\item{even_flag}{a logical value used internally}
\item{nSpecies}{the number of species in the design}

The \code{DI_data_E_AV}, \code{DI_data_FG} and \code{DI_data_ADD} functions return a named list including one or more of the components below (depending on the data manipulation function):
\item{AV}{the AV variable}
\item{E}{the E variable}
\item{ADD_vars}{the variables used in the ADD model}
\item{ADD_theta}{the variables used in the ADD model when theta is estimated}
\item{FG}{the variables used in the FG model}
\item{even_flag}{a logical value used internally}
\item{P_int_flag}{a logical value used internally}

The \code{DI_data_fullpairwise} function returns a matrix with the pairwise interactions between species.
}
\references{
Connolly J, T Bell, T Bolger, C Brophy, T Carnus, JA Finn, L Kirwan, F Isbell, J Levine, A \enc{Lüscher}{}, V Picasso, C Roscher, MT Sebastia, M Suter and A Weigelt (2013) An improved model to predict the effects of changing biodiversity levels on ecosystem function. Journal of Ecology, 101, 344-355. 
 
Kirwan L, A \enc{Lüscher}{}, MT Sebastia, JA Finn, RP Collins, C Porqueddu, A Helgadottir, OH Baadshaug, C Brophy, C Coran, S Dalmannsdottir, I Delgado, A Elgersma, M Fothergill, BE Frankow-Lindberg, P Golinski, P Grieu, AM Gustavsson, M \enc{Höglind}{}, O Huguenin-Elie, C Iliadis, M \enc{Jørgensen}{}, Z Kadziuliene, T Karyotis, T Lunnan, M Malengier, S Maltoni, V Meyer, D Nyfeler, P Nykanen-Kurki, J Parente, HJ Smit, U Thumm, & J Connolly (2007) Evenness drives consistent diversity effects in intensive grassland systems across 28 European sites. Journal of Ecology, 95, 530-539.  
 
Kirwan L, J Connolly, JA Finn, C Brophy, A \enc{Lüscher}{}, D Nyfeler and MT Sebastia (2009) Diversity-interaction modelling - estimating contributions of species identities and interactions to ecosystem function. Ecology, 90, 2032-2038. 
}
\author{
Rafael A. Moral, John Connolly and Caroline Brophy
} 
%\note{
%further notes
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DI}}
\code{\link{autoDI}}

Other examples using the data manipulation functions: 
The \code{\link{Bell}} dataset examples. 
The \code{\link{sim2}} dataset examples. 
The \code{\link{sim3}} dataset examples. 
The \code{\link{sim4}} dataset examples. 
The \code{\link{sim5}} dataset examples. 
The \code{\link{Switzerland}} dataset examples. 
}
\examples{
################################
  
#### Data manipulation for the Switzerland dataset
  
## Load the Switzerland data
  data(Switzerland)
  
## Check that the proportions sum to 1 (required for DI models)
## p1 to p4 are in the 4th to 7th columns in Switzerland
  Switzerlandsums <- rowSums(Switzerland[4:7])
  summary(Switzerlandsums)
  
  
## Create new interaction variables and incorporate them into a new data frame Switzerland2.
## Switzerland2 will contain the new variables:  AV, E, p1_add, p2_add, p3_add, p4_add, 
##  bfg_G_L, wfg_G and wfg_L.
  newlist <- DI_data_prepare(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                             data = Switzerland)
  Switzerland2 <- data.frame(newlist$newdata, newlist$FG)
  
  
## Create new interaction variables and incorporate them into a new data frame Switzerland3.
## Use theta = 0.5.
  newlist <- DI_data_prepare(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                             data = Switzerland, theta = 0.5)
  Switzerland3 <- data.frame(newlist$newdata, newlist$FG)
## Add "_theta" to the new interaction variables to differentiate from when theta = 1
  names(Switzerland3)[12:20] <- paste0(names(Switzerland3)[12:20], "_theta") 
  
  
#### The various interactions can also be added to a new dataset individually:
  
## Create the average pairwise interaction and evenness variables
##  and store them in a new data frame called Switzerland4.
## Switzerland4 will contain the new variables: AV, E
  newlist <- DI_data_E_AV(prop = c("p1","p2","p3","p4"), data = Switzerland)
  Switzerland4 <- data.frame(Switzerland, "AV" = newlist$AV, "E" = newlist$E)
  
## Create the functional group variables and add them to Switzerland4.
## In the FG names vector: G stands for grass, L stands for legume.
## Switzerland4 will contain: bfg_G_L, wfg_G and wfg_L
  newlist <- DI_data_FG(prop = 4:7, FG = c("G","G","L","L"), data = Switzerland)
  Switzerland4 <- data.frame(Switzerland4, newlist$FG)
  
## Create the additive species variables and add them to Switzerland4.
## Switzerland4 will contain the new variables: p1_add, p2_add, p3_add and p4_add.
  newlist <- DI_data_ADD(prop = c("p1","p2","p3","p4"), data = Switzerland)
  Switzerland4 <- data.frame(Switzerland4, newlist$ADD)
  
## Create all pairwise interaction variables and add them to Switzerland4.
## Switzerland5 will contain the new variables: p1.p2, p1.p3, p1.p4, p2.p3, p2.p4, p3.p4.
  newlist <- DI_data_fullpairwise(prop = c("p1","p2","p3","p4"), data = Switzerland)
  Switzerland4 <- data.frame(Switzerland4, newlist)
  
################################
  

################################ 
  
#### Short worked example (as illustrated the Details section)
  
## Create a dataframe
  p1 <- c(0.1, 0.25)
  p2 <- c(0.2, 0.25)
  p3 <- c(0.3, 0.25)
  p4 <- c(0.4, 0.25)
  minidataset1 <- data.frame(p1,p2,p3,p4)
  
## Check the rows sum to 1
  rowSums(minidataset1[1:4]) 
  
## Create the interaction variables, assume two functional groups and theta = 1
  newlist <- DI_data_prepare(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"),
                             data = minidataset1)
  minidataset2 <- data.frame(newlist$newdata, newlist$FG)
  
## Create the interaction variables, assume two functional groups and theta = 0.5
  newlist <- DI_data_prepare(prop = c("p1","p2","p3","p4"), FG = c("G","G","L","L"), 
                             y = "response", data = minidataset1, theta = 0.5)
  minidataset3 <- data.frame(newlist$newdata, newlist$FG, check.names = FALSE)
## Add "_theta" to the new interaction variables to differentiate from when theta = 1
  names(minidataset3)[8:16] <- paste0(names(minidataset3)[8:16], "_theta")   
  
################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
