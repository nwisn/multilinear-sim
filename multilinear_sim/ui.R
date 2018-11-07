#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("superhero"),

  # Application title
  titlePanel("Monte-Carlo Comparison of Regression Models"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
       p("Generate data Y from a high-dimensional linear model, with features X sampled from a skew-normal distribution, and random error from a normal distribution. Then randomly discard information (features) and build predictive models using OLS and shrinkage estimators, with and without an inverse-density weighting scheme based on the distribution of Y."),
       numericInput("p",
                   "Dimensionality of generative model:",
                   min = 1,
                   max = 5000,
                   step = 1,
                   value = 200),
       numericInput("frac",
                    "Fraction of features observed:",
                    min = 0,
                    max = 1,
                    step = 0.1,
                    value = 0.6),
       numericInput("epsilon",
                    "Random error:",
                    min = 0,
                    max = 100,
                    step = .1,
                    value = 12),
       numericInput("alpha",
                    "Alpha (skew-normal):",
                    min = 0,
                    max = 100,
                    step = .1,
                    value = 0),
       numericInput("n",
                    "Sample size:",
                    min = 10,
                    max = 10000,
                    step = 100,
                    value = 500),
       numericInput("randomseed",
                    "Random seed:",
                    min = 1,
                    max = 1000000000,
                    step = 1,
                    value = 17)

    ),

    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("figure", width = 600, height = 900)
    )
  )
))
