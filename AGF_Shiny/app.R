
library(shiny)
library(tidyverse)

#### Reading the data files
wideFile_Fitness <- read.csv("/Users/whitlock/Desktop/AGF_Shiny/supersetMeansRelFitnessA06.csv", stringsAsFactors = FALSE)
wideFile_Replacement <- read.csv("/Users/whitlock/Desktop/AGF_Shiny/supersetMeansReplacementA06.csv", stringsAsFactors = FALSE)


### Converting the wide file into long format  --NOTE: It may be worth embedding
### this into the reactive part (after filtering) to speed up the start-up.
# Fitness df
gen =  gather(wideFile_Fitness, key=genChar, value = relFitness, "g001":"g100", factor_key=TRUE)
gen$genChar=substr(gen$genChar,2,4)
gen$Generation=as.numeric(gen$genChar)

#Genome replacement df
genRep =  gather(wideFile_Replacement, key=genChar, value = hybridIndex, "g001":"g100", factor_key=TRUE)
genRep$genChar=substr(genRep$genChar,2,4)
genRep$Generation=as.numeric(genRep$genChar)


### User interface

ui <- fluidPage(
    titlePanel(HTML("<h2>Assisted Gene Flow (AGF)</h2><h3>Fitness and proportion local genotypes remaining </h3>")),
    sidebarLayout(
        sidebarPanel(
            radioButtons("responseVar", "Response variable",
                         choices = c( "Fitness", "Local genes remaining"),
                         selected = "Fitness"),
            radioButtons("popSize", "Population size, N",
                         choices = c( 1000, 10000),
                         selected = 10000,
                         inline=TRUE),
            radioButtons("nbLoci_PAC", "Number of pre-adapted loci",
                         choices = c(0, 1, 5, 50),
                         selected = 5,
                         inline=TRUE),
            radioButtons("nbLoci_MAC", "Number of locally maladaptive loci",
                         choices = c(0, 1, 5, 50),
                         selected = 0,
                         inline=TRUE),
            radioButtons("nbLoci_OD", "Number of pairs of outbreeding depression loci",
                         choices = c(2, 10, 100),
                         selected = 10,
                         inline=TRUE),
            # radioButtons("numPulses", "Number of pulses of introduction",
            #              choices = c(1,5),
            #              selected = 1,
            #              inline=TRUE),
           radioButtons("pulseInterval", "Introduction pulse frequency", 
                             choiceNames = list("Every generation", "Every 2 gens","Every 4 gens"),
                             choiceValues = list(1,2,4),
                             selected = 1),
                
        
            radioButtons("migrationRate", HTML("<p>T<sub>f</sub>, proportion introduced</p>"),
                         choices = c(0.005, 0.05, .5),
                         selected = 0.05,
                         inline=TRUE),
            radioButtons("fit_PAC_param",
                         HTML("<p >&Delta;<sub>PA</sub>, Maximum fitness change, preadapted loci</p>"),
                         choices = c(0.1, 0.5),
                         selected = 0.1,
                         inline=TRUE),
            conditionalPanel(
                "input.nbLoci_MAC>0",
                radioButtons("fit_MAC_param",
                             HTML("<p >&Delta;<sub>MA</sub>, Maximum fitness change, maladapted loci</p>"),
                         choices = c(0.1, 1, 10),      
                         inline=TRUE),
                selected = 0.1),
            
            sliderInput("generations", "Generations", min = 0, max = 100,
                    value = c(0, 40) ),
            uiOutput("GeneticRescue")
        ),
        
       mainPanel( HTML("<p><b>Single introduction</b></p>"),
                  plotOutput("plotOnePulse"),
                  HTML("<p><b>Five pulses of introduction</b></p>"),
                  plotOutput("plotFivePulse")
                  ),
    ),
####  Footer note    
    HTML("<div id = \"credits\">Developed as supplementary material for  <a href=\"https://www.biorxiv.org/\">a paper exploring the fitness effects of AGF</a>  
        by Jared Grummer, Tom Booker, Rémi Matthey-Doret, Pirmin Nietlisbach, Andrea Thomaz, and Mike Whitlock at the University of British Columbia. </p></div>")
)

server <- function(input, output) {
    response = reactive({input$responseVar})
    
### filtering the data set to match the input choices
    filteredFitness <- reactive({
        #pI=if(input$numPulses==1) 0 else input$pulseInterval
        singleSMAC=if(input$nbLoci_MAC==0) 0 else input$fit_MAC_param
        singleSPAC=if(input$nbLoci_PAC==0) 0 else input$fit_PAC_param
        gen %>%
            filter(Generation >= input$generations[1],
                   Generation <= input$generations[2],
                   N == input$popSize,
                   #numPulses == input$numPulses | numPulses == 1,
                   #pulseInterval==input$pulseInterval,
                   totalReplacementFraction==input$migrationRate,
                   singleSToRescale_PAC == singleSPAC,
                   singleSToRescale_MAC == singleSMAC,
                   nbLoci_MAC == input$nbLoci_MAC,
                   nbLoci_PAC == input$nbLoci_PAC,
                   nbPairsLoci_OD == input$nbLoci_OD | nbPairsLoci_OD ==0,
                   #dominanceCoef_MAC==0.5
            )
    })
    
    filteredReplacement <- reactive({
       # pI=if(input$numPulses==1) 0 else input$pulseInterval
        singleSMAC=if(input$nbLoci_MAC==0) 0 else input$fit_MAC_param
        singleSPAC=if(input$nbLoci_PAC==0) 0 else input$fit_PAC_param
        genRep %>%
            filter(Generation >= input$generations[1],
                   Generation <= input$generations[2],
                   N == input$popSize,
                  #numPulses == input$numPulses,
                   #pulseInterval==pI,
                   totalReplacementFraction==input$migrationRate,
                   singleSToRescale_PAC == singleSPAC,
                   singleSToRescale_MAC == singleSMAC,
                   nbLoci_MAC == input$nbLoci_MAC,
                   nbLoci_PAC == input$nbLoci_PAC,
                   nbPairsLoci_OD == input$nbLoci_OD | nbPairsLoci_OD ==0,
                   #dominanceCoef_MAC==0.5
            )
    })
    
    minFit = reactive({min(filteredFitness()$relFitness)})
    maxFit = reactive({max(filteredFitness()$relFitness)})
    minRep = reactive({min(filteredReplacement()$hybridIndex)})
    maxRep = reactive({max(filteredReplacement()$hybridIndex)})

####  Making the plot    

    
    output$plotOnePulse <- renderPlot({
        if (is.null(filteredFitness())) {
            return()
        }
        if(input$responseVar=="Fitness"){
            ggplot( filteredFitness() %>% filter(numPulses==1), aes(x = Generation, y = relFitness, col = DeltaOD)) +
                geom_line(size = 1.5)+
                geom_hline(yintercept=1, linetype="longdash", color = "grey", size = 0.5) +
                labs(x="Generations", y = "Mean relative fitness of recipient population",
                     colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
                scale_color_manual(values=c("black", "#2c7fb8", "#e6550d")) +
                ylim(minFit(),maxFit())+
                theme_classic()+
                theme(text = element_text(size=16))+
                theme(aspect.ratio=2/3)
            
        } else
            if(input$responseVar=="Local genes remaining"){
                ggplot( filteredReplacement()%>% filter(numPulses==1), aes(x = Generation, y = hybridIndex, col = DeltaOD)) +
                    geom_line(size = 1.5)+
                    
                    labs(x="Generations", y = "Proportion resident genotype \n remaining in recipient population",
                         colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
                    scale_color_manual(values=c("black", "#2c7fb8", "#e6550d")) +
                    ylim(minRep(),maxRep())+
                    theme_classic()+
                    theme(text = element_text(size=16))+
                    theme(aspect.ratio=2/3)
            }
    })
    
    output$plotFivePulse <- renderPlot({
        if (is.null(filteredFitness())) {
            return()
        }
        if(input$responseVar=="Fitness"){
            ggplot( filteredFitness()%>% filter(numPulses==5) %>% filter(pulseInterval==input$pulseInterval), 
                    aes(x = Generation, y = relFitness, col = DeltaOD)) +
                geom_line(size = 1.5)+
                geom_hline(yintercept=1, linetype="longdash", color = "grey", size = 0.5) +
                labs(x="Generations", y = "Mean relative fitness of recipient population",
                     colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
                scale_color_manual(values=c("black", "#2c7fb8", "#e6550d")) +
                ylim(minFit(),maxFit())+
                theme_classic()+
                theme(text = element_text(size=16))+
                theme(aspect.ratio=2/3)
        } else 
            if(input$responseVar=="Local genes remaining"){
                ggplot( filteredReplacement()%>% filter(numPulses==5) %>% filter(pulseInterval==input$pulseInterval), 
                        aes(x = Generation, y = hybridIndex, col = DeltaOD)) +
                    geom_line(size = 1.5)+
                    
                    labs(x="Generations", y = "Proportion resident genotype \n remaining in recipient population",
                         colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
                    scale_color_manual(values=c("black", "#2c7fb8", "#e6550d")) +
                    ylim(minRep(),maxRep())+
                    theme_classic()+
                    theme(text = element_text(size=16))+
                    theme(aspect.ratio=2/3)
            }
        
    })
    
    # output$plotReplacement <- renderPlot({
    #     if (is.null(filteredReplacement())) {
    #         return()
    #     }
    #     ggplot( filteredReplacement(), aes(x = Generation, y = hybridIndex, col = DeltaOD)) +
    #         geom_line(size = 1.5)+
    # 
    #         labs(x="Generations", y = "Proportion genotype remaining in recipient population",
    #              colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
    #         scale_color_manual(values=c("black", "#2c7fb8", "#e6550d")) +
    #         ylim(minRep(),maxRep())+
    #         theme_classic()+
    #         theme(text = element_text(size=16))
    # 
    # })
}


shinyApp(ui = ui, server = server)
