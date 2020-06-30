
## Input:
# p: the vector pf probability that the player 1, 2, ..., n wins in a single trial
# k1: the vector regarding the number of trials that player 1, 2, ..., n have won when interrupted
# t: the number of successes of trials required to win the game
# simulation: the number of games conducted to get the simulated probability that player 1 wins the game

## Output:
# The simulated and calculated value of the probability that player 1 wins. And the relative different between them.
library(MGLM)
library(HDInterval)
library(DT)
library(MCMCpack)
#source("combinations.R")



compare_function1 <- function(p,m1,n1,t,simulation){
  
  # first estimating the probability that player 1 wins by simulating the game for a given number of times.
  m <- t - m1 # m: the number of success that player 1 should have to win the game
  n <- t - n1 # n: the number of success that player 2 should have to win the game
  
  record_game <- vector()
  for (i in 1:simulation){
    # copy m and n
    record_trial <- rbind(m,n)
    # keep running the game until one of them wins
    while (record_trial[1,]*record_trial[2,]!= 0){
      record_trial <- record_trial - rmultinom(1,1,c(p,1-p)) #updating the recorded result 
    }
    
    # if record_trial[1,] = 0, record 1,player 1 wins this single trial; 
    # if record_trial[2,] = 0, record 0, player 2 wins 
    record_game[i] <- record_trial[2,]^record_trial[1,]
    
  }
  p_simulated <- sum(record_game)/simulation    # estimated probability that player 1 wins
  
  # plot the simulated probability over the number of times the simulation is performed
  p_cumulative <- record_game
  for (i in 2:simulation){
    p_cumulative[i] <- (p_cumulative[i] + p_cumulative[i-1]*(i-1))/i
  }
  #plot(p_cumulative,
  #main="Simulated probability of player 1 winning over simulations", 
  #xlab = "Number of simulations",ylab="Probability of winning",type="l",
  #xlim = c(0, simulation), ylim = c(0, 1),)
  
  # calculating the probabiltiy that player 1 wins from negative binomial distribution
  p_calculated <- pnbinom(n-1,m,p)  
  
  # difference between simulated and calculated probability
  dif = abs(p_simulated - p_calculated) 
  return (list(p_simulated,p_calculated,dif,p_cumulative))
}

compare_function2 <- function(k1, t, p, simulation){
  # first estimating the probability that player 1 wins
  k <- t-k1 # k: the vector regarding the number of more successes that player 1, 2, ..., n should have to win the game
  
  # normalize the probability if they do not sum to one
  p <- p/sum(p)
  
  record_game <- vector()# initializing the vector to record the results of all the games
  for (i in 1:simulation){
    # copy the trials players need to succeed and update after each trial
    record_trial <- k # each row stands for simulations
    # keep running the game until one record of the trial reaches 0
    while (prod(record_trial) != 0){
      record_trial <- record_trial - t(rmultinom(1,1,p)) #updating the recorded result
    }
    # transfer the record of trial to the record of game
    record_game[i] <- which(record_trial == 0, arr.ind = TRUE)[1,"col"]
  }
  
  # estimated probability that each player wins
  p_simulated <- vector()
  for (i in 1:length(p)){
    p_simulated[i] <- length(which(record_game == i))/simulation 
  }
  
  p_cumulative <- matrix(0, nrow = length(p), ncol = simulation)
  p_cumulative[record_game[1],1] <- 1
  for (i in 1:length(p)){    
    for (j in 2:simulation){
      
      if (record_game[j] == i){
        p_cumulative[i, j] <- (p_cumulative[i, j-1]*(j-1)+1)/j
      }else{
        p_cumulative[i, j] <- p_cumulative[i, j-1]*(j-1)/j
      }
    }
  }

  ## Calculating the probability using negative multinomial distribution and resursive function
  p_calculated <- 0
  for (l in 1:(prod(k)/k[1])){
    p_calculated <- p_calculated + exp(dnegmn(Y = combinations(k)[l,2:length(k)],
                                                beta = combinations(k)[1], prob = p[2:length(p)]))
  }
  
  # calculating difference between the simulated and the calculated probability of player 1 winning
  dif <- abs(p_simulated - p_calculated)
  return (list(p_simulated, p_calculated, dif, p_cumulative))
}


#l <- compare_function2(c(1,1,1), 2, c(0.2,0.4,0.4), 10)
#l

combinations <- function(k){
  if (length(k) == 3){ # the number of players
    #k <- c(1,2,3)
    num <- prod(k)/k[1] #total number of sets
    #sets <- vector(mode = "list", length = num)
    sets <- matrix(nrow = num, ncol = length(k))
    for (i in 0:(k[2]-1)){
      for (j in 0:(k[3]-1)){
        sets[i*k[3]+j+1,] <- c(k[1],i,j)
      }
    }
    return(sets)
  }else{
    #k <- c(1,2,3,4)
    #x1 <- sets
    x1 <- combinations(k[1:(length(k) - 1)]) # k without k[length(k)]
    x2 <- k[length(k)]  
    y <- matrix(ncol = length(k), nrow = prod(k)/k[1])
    for (i in 1:length(x1[,1])){
      for (j in 1:k[length(k)]){
        y[(i-1)*(k[length(k)])+j,] <- append(x1[i,],j-1)
      }
    }
    return(y)
  }
} 




library(shiny)

# Create the shinyapp
ui <- fluidPage(
  titlePanel("Game of Chance for N Players"),
  sidebarPanel(width = 3,
               numericInput(inputId = "n", label = "The number of players",value = 3), 
               textInput(inputId = "k", label = "The number of points for each player when interrupted (comma delimited)", 
                         value = "1,1,1"),
               numericInput(inputId = "t", label = "The number of points required to win the game",value = 2), 
               textInput(inputId = "p", label = "For every play, the probability that each player wins the point (comma delimited)", 
                         value = "0.3,0.2,0.5"),
               numericInput(inputId = "s", label = "The number of simulated games", 
                            value = 1000),
               checkboxInput(inputId = "check", label = "95% credible interval (highest posterior density)") 
  ),
  mainPanel(
    # textOutput("vector"),
    textOutput("warning"),
    DTOutput('table'),
    plotOutput("simulation"),
    #tags$style(type="text/css",
    #           ".shiny-output-error { visibility: hidden; }",
    #           ".shiny-output-error:before { visibility: hidden; }"
    #)
    
  )
)




server <- function(input, output) {
  
  
  ## user input: vector of the number of points won
  k <-  reactive({
    as.numeric(unlist(strsplit(input$k,",")))
  })
  
  ## user input: vector of probability of winning each point
  p <- reactive({
    as.numeric(unlist(strsplit(input$p,",")))/sum(as.numeric(unlist(strsplit(input$p,","))))
  })
  
  ## output of compare_fonction1, when there are two players
  result1 <- reactive({
    compare_function1(p()[1],k()[1],k()[2],input$t,input$s)
  })
  
  ## output of compare_function2, when there are three or more players
  result2 <- reactive({
    compare_function2(k(), input$t, p(), input$s)
  })
  #try(if(max(k()) < input$t) stop("Warning:"))
  # warning output if any player has won the game
  output$warning <- renderText({ 
    
    validate(
      need(max(k()) < input$t, paste("Warning: Player", which(k() == max(k())),
                                     " has already won the game. Adjust the inputs!")),
      
      #need(max(k()) < input$t, stop("Warning: Player", which(k() == max(k())),
      #                               " has already won the game. Adjust the inputs!")),
      
      need(sum(c(k(),p()) > 0) == length(c(k(),p())),"Warning: No negative input values!"),
      need(input$n >= 2,"Warning: The number of players must be at least 2!"),
      need(sum(as.numeric(unlist(strsplit(input$p,",")))) == 1, 
           "Warning: The probability does not sum to 1. Already normalized!")
      
    )  
  })
  
  
 
  
  output$table <- renderDT({
    # adapt the header
    sketch <- 
      htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(rowspan = 2, 'Players'),
            th(rowspan = 2, 'Pr(win a point)'),
            th(colspan = 2, 'Pr(win the game)')
          ),
          tr(
            lapply(c('Analytical', 'Simulated'), th)
          )
        )
      ))
    
    if (input$n == 2 & max(k()) < input$t){
      # specify the output from compara_fuction1 in a data frame 
      data <- data.frame("Player" = 1:2, 
                         "Prob(winning a trial)" = p(),
                         "Analytical Prob of winning the game" = c(result1()[[2]],1-result1()[[2]]),
                         "Simulated Prob of winning the game" = c(result1()[[1]],1-result1()[[1]]),
                         check.names = FALSE) 
      datatable(data, container = sketch, rownames = FALSE, caption = "Summary Table", 
                options = list(dom = 't'))
    }else if (input$n >= 3 & max(k()) < input$t){
      ## a vector of analytical p for all players, calculated by switching with player 1
      Analytical_Prob <- vector()
      for (i in 1:length(p())){
        k_copy <- replace(k(), c(1, i), k()[c(i, 1)])
        p_copy <- replace(p(), c(1, i), p()[c(i, 1)])
        Analytical_Prob[i] <- compare_function2(k_copy, input$t, p_copy, input$s)[[2]] 
      }
      # specify the output from compara_fuction1 in a data frame 
      data <- data.frame("Player" = 1:length(p()), 
                         "Prob(winning a trial)" = p(),
                         "Analytical Prob of winning the game" = Analytical_Prob,
                         "Simulated Prob of winning the game" = result2()[[1]],
                         check.names = FALSE) 
      datatable(data, container = sketch, rownames = FALSE, caption = "Summary Table", 
                options = list(dom = 't'))
    }
    
  })
  
  output$simulation <- renderPlot({
    if (input$n == 2 & max(k()) < input$t){
      plot(x <- c(1:input$s), y = rep(result1()[[2]], input$s),
           main="Probability of Player 1 Winning", cex.main = 2,
           xlab = "Number of Simulated Games",ylab="Pr(Winning the Game)",cex.lab = 1.5,type="l",
           xlim = c(0, input$s), ylim = c(0, 1), lwd = 2, col = "red") 
      
      lines(x = c(0.73*input$s,0.79*input$s), y = c(0.9,0.9), pch = 21, lwd = 2, cex = 1.5)
      text(0.8*input$s, 0.9, "Simulated probability", cex = 1.2, font = 1, adj = 0)
      lines(x = c(0.73*input$s,0.79*input$s), y = c(0.8,0.8), pch = 21, lwd = 2, cex = 1.5,col = "red")
      text(0.8*input$s, 0.8, "Analytical probability", cex = 1.2, font = 1, adj = 0,col = "red")
      
      ##Testing
      #input <- data.frame("s" = 6)
      #SimulResult <- c(1,0,1,0,1,0)
      
      ## Credibility interval (highest posterior density interval)
      SimulResult <- result1()[[4]] # store the simulated result
      SimulMatrix <- matrix(0, nrow = 1000, ncol = input$s) # the matrix of samples from posterior distribution based on simulated result
      
      for (i in 1:input$s){  
        SimulMatrix[, i] <- rbeta(1000, SimulResult[i]*i+1, i-SimulResult[i]*i+1)
      }
      CredInt <- apply(SimulMatrix, 2, hdi) # record the credibility interval
      
      y.upper <- CredInt[1,]
      y.lower <- CredInt[2,]
      #lines(y.upper)
      #lines(y.lower)
      if (input$check){
        polygon(c(1:input$s,input$s:1), c(y.upper, rev(y.lower)), col = "lightsteelblue", border = NA)
        lines(x = c(1:input$s), y = rep(result1()[[2]], input$s), col = "red")
      }
      ## simulated probability
      lines(x <- c(1:input$s), y = result1()[[4]]) 
    }else if(input$n >= 3 & max(k()) < input$t){
      par(las = 1)
      plot(result2()[[4]][1,],
           main="Probability of Player 1 Winning", cex.main = 2,
           xlab = "Number of Simulated Games",ylab="Pr(Winning the Game)",
           cex.lab = 1.5,type="l",
           xlim = c(0, input$s), ylim = c(0, 1))
      lines(x = c(1:input$s), y = rep(result2()[[2]], input$s),col = "red")
      lines(x = c(0.73*input$s,0.79*input$s), y = c(0.9,0.9), pch = 21, lwd = 2, cex = 1.5)
      text(0.80*input$s, 0.9, "Simulated probability", cex = 1.2, font = 1, adj = 0)
      lines(x = c(0.73*input$s,0.79*input$s), y = c(0.8,0.8), pch = 21, lwd = 2, cex = 1.5,col = "red")
      text(0.80*input$s, 0.8, "Analytical probability", cex = 1.2, font = 1, adj = 0,col = "red")
      
      
      ## Credibility interval (highest posterior density interval)
      SimulResult <- result2()[[4]][1,] # store the simulated result
      SimulMatrix <- matrix(0, nrow = 1000, ncol = input$s) # the matrix of samples from posterior distribution based on simulated result
      
      
      
      for (i in 1:input$s){  
        #SimulMatrix[, i] <- rbeta(1000, SimulResult[i]*i+1, i-SimulResult[i]*i+length(p())-1)
        SimulMatrix[, i] <- rdirichlet(1000, result2()[[4]][,i]*i+1)[,1]
      }
      CredInt <- apply(SimulMatrix, 2, hdi) # record the credibility interval
      
      y.upper <- CredInt[1,]
      y.lower <- CredInt[2,]
      #lines(y.upper)
      #lines(y.lower)
      if (input$check){
        polygon(c(1:input$s,input$s:1), c(y.upper, rev(y.lower)), col = "lightsteelblue", border = NA)
        lines(x = c(1:input$s), y = rep(result2()[[2]], input$s), col = "red")
      }
      ## simulated probability
      lines(x <- c(1:input$s), y = result2()[[4]][1,]) 
    }
    
  })
}

shinyApp(ui = ui, server = server)

