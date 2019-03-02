

shinyServer(function(input, output, session){
  InputPathway <- eventReactive(
    input$pathways,
    {
      switch(input$pathways,
             "Apoptosis" = Apoptosis_gene_string,
             "DNA Repair" = DNA_string_genes,
             "Cell Cycle" = CellCycle_gene_string,
             "Hedgehog signaling" = hedgehog_gene_string,
             "TCA Cycle" = TCA_gene_string,
             "Adherens Junction" = Adherens_gene_string,
             "Pancreatic Cancer" = Pancreatic_gene_string,
             "ABC Transporters" = ABC_gene_string)
    }
  )
  
  observeEvent(
    {
      input$Tables
      input$Regulation
      input$pathways
    },
    {
      reac_table_1 <- reactive(table0 <- deSets[[paste(input$Tables, input$Regulation, sep=".")]])
      pathway <- InputPathway()
      
      genes.string <- str_c(unlist(reac_table_1()$gene), sep = ",")
      intection.pathway <- intersect(pathway, genes.string)
      
      
      if(length(intection.pathway) != 0){
        positions <- match(pathway, genes.string)
        #positions <- na.omit(positions)
        output$table1 <- renderDataTable({datatable(na.omit(reac_table_1()[positions, ]))})
        
        output$Download.tab5.1 <- downloadHandler(
          filename = function(){
            paste(input$pathways,"-",input$Tables,"-",input$Regulation,".csv" ,sep="")
          },
          content = function(file){
            write.csv(na.omit(reac_table_1()[positions, ]), file, col.names = FALSE, row.names = FALSE, quote = FALSE
            )
          },
          contentType = "csv"
        )
      }
      
    }
  )# observeevent1 ends
  
  
  
  
  
  
  
  
  # For the loading file
  filedata <- eventReactive(
    
    input$Upload.file,
    {
      infile <- input$Upload.file
      if(is.null(infile)){return(NULL)
      }else{
        file <-  read.csv(infile$datapath, header = FALSE)
        datagene <- as.character(file[ ,1])
      }
    }
  )
  
  
  observeEvent(
    {
      input$Upload.file
      input$Tables
      input$Regulation
    },
    {
      file.genes <- filedata()
      dataset.genes <-str_c(unlist(deSets[[paste(input$Tables, input$Regulation, sep=".")]][ ,1]), sep = ",")
      exist.genes <- intersect(file.genes, dataset.genes)
      
      
      genes.position <- match(file.genes, dataset.genes)
      output$gene.text <- renderDataTable({
        datatable(na.omit(deSets[[paste(input$Tables, input$Regulation, sep=".")]][genes.position, ]))
      })
      output$Download.tab5.2 <- downloadHandler(
        filename = function(){
          paste("genelist","-",input$Tables,"-",input$Regulation,".csv" ,sep="")
        },
        content = function(file){
          write.csv(na.omit(deSets[[paste(input$Tables, input$Regulation, sep=".")]][genes.position, ]),
                    file, col.names = FALSE, row.names = FALSE, quote = FALSE)
        },
        contentType = "csv"
      )
      
      
    }
  )
  
  
  
  
  
  observeEvent(
    input$check_sort_1,
    {
      
      newTable1 <- reactive(Table1 <- deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1, ])
      
      
      output$topDEgenes1 <- renderDataTable({
        datatable(newTable1())
      })
      
      output$downloadData1 <- downloadHandler(
        
        
        filename = function(){
          paste(input$deTest1,"-",input$direction1,".csv", sep="")
          
        },
        
        content = function(file){
          
          write.csv(newTable1(), file, col.names = FALSE, row.names = FALSE, quote = FALSE)
          
        },
        contentType = "csv"
      )
      
    }
  )# first observe
  
  observeEvent(
    input$check_sort_2,
    {
      newTable2 <- reactive(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2, ])
      
      output$topDEgenes2 <- renderDataTable({
        datatable(newTable2())
      })
      
      output$downloadData2 <- downloadHandler(
        filename = function(){
          paste(input$deTest2,"-", input$direction2,".csv", sep="")
          
        },
        
        content = function(file){
          
          write.csv(newTable2(), file, col.names = FALSE, row.names = FALSE, quote = FALSE)
          
        }
        
        
      )
      
    }
  )
  
  observeEvent(
    input$check_sort_3,
    {
      
      newTable3 <- reactive(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1, ])
      output$topDEgenes3 <- renderDataTable({
        datatable(newTable3())
      })
      
      output$downloadData1 <- downloadHandler(
        
        filename = function(){
          paste(input$deTest1,"-", input$direction1,"csv", sep="")
          
        },
        
        content = function(file){
          
          write.csv(newTable3(), file, col.names = FALSE, row.names = FALSE, quote = FALSE )
          
        }
      )
      
    }
  )
  
  observeEvent(
    input$check_sort_4,
    {
      newTable4 <- reactive(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2, ])
      output$topDEgenes4 <- renderDataTable({
        datatable(newTable4())
      })
      
      output$downloadData2 <- downloadHandler(
        filename = function(){
          paste(input$deTest2,"-", input$direction2,".csv", sep="")
          
        },
        
        content = function(file){
          
          write.table(newTable4(), file, col.names = FALSE, row.names = FALSE, quote = FALSE )
        }
      )
      
    }
  )
  
  
  
  
  # Draw diagram
  
  output$vennPlot <- renderPlot({
    other.rec1 <- reactiveValues(reac1 = input$deTest1)
    other.rec2 <- reactiveValues(reac2 = input$deTest2)
    other.rec3 <- reactiveValues(reac3 = input$direction1)
    other.rec4 <- reactiveValues(reac4 = input$direction2)
    other.rec5 <- reactiveValues(reac5 = input$number1)
    other.rec6 <- reactiveValues(reac6 = input$number2)
    
    
    # Using if else statement to eliminate some coincidences
    if( other.rec1$reac1 == other.rec2$reac2  & other.rec3$reac3 == other.rec4$reac4 & other.rec5$reac5 == other.rec6$reac6){
      return()
      #print("they are same!")
      
    }else{
      
      # Reactive the areas in venn diagram
      Area1 <- reactiveValues(genes1 = as.character(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1,1]))
      Area2 <- reactiveValues(genes2 = as.character(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2,1]))
      # area1, area2, cross.area
      lengthArea1 <- reactiveValues(length1 = length(Area1$genes1))
      lengthArea2 <- reactiveValues(length2 = length(Area2$genes2))
      overlap.area <- reactiveValues(length3 = length(intersect(Area1$genes1 , Area2$genes2)))
      
      draw.pairwise.venn(
        
        area1 = lengthArea1$length1,
        area2 = lengthArea2$length2,
        cross.area = overlap.area$length3,
        category = c(input$deTest1, input$deTest2),
        lty = rep("blank"),
        fill = c("light blue", "pink"),
        
        cex = 2,
        cat.cex = 2,
        cat.pos = c(285, 105),
        cat.dist = 0.09,
        cat.just = list(c(-1, -1), c(1, 1)),
        ext.pos = 30,
        ext.dist = -0.05,
        ext.length = 0.85,
        ext.line.lwd = 2,
        ext.line.lty = "dashed"
        
      )
    }
  })
  
  
  
  
  
  # Download data from venn diagram
  cover_Area <- function(Area1, Area2){
    
    Area1 <- reactiveValues(genes1 = as.character(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1,1]))
    Area2 <- reactiveValues(genes2 = as.character(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2,1]))
    coverArea <- intersect(Area1$genes1 , Area2$genes2)
    return(coverArea)
  }
  
  output$downloadData3 <- downloadHandler(
    filename = function(){
      paste(c(input$deTest1, input$deTest2), c(input$direction1, input$direction2), ".txt", sep="")
      
    },
    
    content = function(file){
      
      write.table(cover_Area(), file, row.names = FALSE, col.names = FALSE, quote = FALSE) 
      
    }
    
  )
  
  
  
  
  
  
  
  
  
  observeEvent(
    input$Button,
    
    {
      Area1 <- reactiveValues(genes1 = as.character(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1,1]))
      Area2 <- reactiveValues(genes2 = as.character(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2,1]))
      overlapGenes <- reactiveValues(overlap = intersect(Area1$genes1 , Area2$genes2))
      ann <- row.names(expressionTable)
      matchgenes <- ann %in% overlapGenes$overlap
      OverlapTable <- expressionTable[matchgenes, ]
      OverlapTable.genes <- row.names(OverlapTable)
      ############################################### 
      observe({updateSelectInput(session, "check_gene", choices = OverlapTable.genes)})
      
      col_include <- length(input$heatmap.1)
      gene_intersect <- input$check_gene %in% rownames(OverlapTable)
      if(isTRUE(gene_intersect) == FALSE){
        gene_intersect <- FALSE
      }
      gene_pos <- rownames(OverlapTable) %in% input$check_gene
      
      
      
      if(col_include != length(colnames(OverlapTable)) & gene_intersect == TRUE){
        OverlapMatrix <- data.matrix(OverlapTable[!(gene_pos),input$heatmap.1])
        output$GeneHeatmap <- renderD3heatmap({
          d3heatmap(OverlapMatrix, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
      else if(col_include != length(colnames(OverlapTable)) & gene_intersect == FALSE){
        OverlapMatrix <- data.matrix(OverlapTable[ ,input$heatmap.1])
        output$GeneHeatmap <- renderD3heatmap({
          d3heatmap(OverlapMatrix, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }#first else if ends
      else if(col_include == length(colnames(OverlapTable)) & gene_intersect == TRUE){
        OverlapMatrix <- data.matrix(OverlapTable[!(gene_pos), ])
        output$GeneHeatmap <- renderD3heatmap({
          d3heatmap(OverlapMatrix, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }# second else if ends
      else if(col_include == length(colnames(OverlapTable)) & gene_intersect == FALSE){
        OverlapMatrix <- data.matrix(OverlapTable)
        output$GeneHeatmap <- renderD3heatmap({
          d3heatmap(OverlapMatrix, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
    }
  )
  
  observeEvent(
    input$Button2,
    
    {
      
      Area1 <- reactiveValues(genes1 = as.character(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1,1]))
      Area2 <- reactiveValues(genes2 = as.character(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2,1]))
      overlapGenes <- reactiveValues(overlap = intersect(Area1$genes1 , Area2$genes2))
      ann2 <- row.names(pancRIP_expression)
      matchgenes2 <- ann2 %in% overlapGenes$overlap
      OverlapTable2 <- pancRIP_expression[matchgenes2, ]
      OverlapTable2.genes <- row.names(OverlapTable2)
      ###############################################
      
      observe({updateSelectInput(session, "check_gene.2", choices = OverlapTable2.genes)})
      
      col_include1 <- length(input$heatmap.2)
      gene_intersect1 <- input$check_gene.2 %in% rownames(OverlapTable2)
      if(isTRUE(gene_intersect1) == FALSE){
        gene_intersect1 <- FALSE
      }
      gene_pos1 <- rownames(OverlapTable2) %in% input$check_gene.2
      
      if(col_include1 != length(colnames(OverlapTable2)) & gene_intersect1 == TRUE){
        OverlapMatrix1 <- data.matrix(OverlapTable2[!(gene_pos1),input$heatmap.2])
        output$GeneHeatmap2 <- renderD3heatmap({
          d3heatmap(OverlapMatrix1, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
      
      else if(col_include1 != length(colnames(OverlapTable2)) & gene_intersect1 == FALSE){
        OverlapMatrix1 <- data.matrix(OverlapTable2[ ,input$heatmap.2])
        output$GeneHeatmap2 <- renderD3heatmap({
          d3heatmap(OverlapMatrix1, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
      
      else if(col_include1 == length(colnames(OverlapTable2)) & gene_intersect1 == TRUE){
        OverlapMatrix1 <- data.matrix(OverlapTable2[!(gene_pos1), ])
        output$GeneHeatmap2 <- renderD3heatmap({
          d3heatmap(OverlapMatrix1, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }# second else if ends
      else if(col_include1 == length(colnames(OverlapTable2)) & gene_intersect1 == FALSE){
        OverlapMatrix1 <- data.matrix(OverlapTable2)
        output$GeneHeatmap2 <- renderD3heatmap({
          d3heatmap(OverlapMatrix1, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
    }
  )
  
  
  observeEvent(
    input$Button3,
    {
      
      Area1 <- reactiveValues(genes1 = as.character(deSets[[paste(input$deTest1, input$direction1, sep=".")]][1:input$number1,1]))
      Area2 <- reactiveValues(genes2 = as.character(deSets[[paste(input$deTest2, input$direction2, sep=".")]][1:input$number2,1]))
      overlapGenes <- reactiveValues(overlap = intersect(Area1$genes1 , Area2$genes2))
      ann3 <- row.names(expression.rnaseq)
      matchgenes3 <- ann3 %in% overlapGenes$overlap
      OverlapTable3 <- expression.rnaseq[matchgenes3, ]
      OverlapTable3.genes <- row.names(OverlapTable3)
      
      ####################################
      observe({updateSelectInput(session, "check_gene.3", choices = OverlapTable3.genes)})
      
      col_include2 <- length(input$heatmap.3)
      gene_intersect2 <- input$check_gene.3 %in% rownames(OverlapTable3)
      if(isTRUE(gene_intersect2) == FALSE){
        gene_intersect2 <- FALSE
      }
      gene_pos2 <- rownames(OverlapTable3) %in% input$check_gene.3
      
      
      
      if(col_include2 != length(colnames(OverlapTable3)) & gene_intersect2 == TRUE){
        OverlapMatrix2 <- data.matrix(OverlapTable3[!(gene_pos2),input$heatmap.3])
        output$GeneHeatmap3 <- renderD3heatmap({
          d3heatmap(OverlapMatrix2, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
      else if(col_include2 != length(colnames(OverlapTable3)) & gene_intersect2 == FALSE){
        OverlapMatrix2 <- data.matrix(OverlapTable3[ ,input$heatmap.3])
        output$GeneHeatmap3 <- renderD3heatmap({
          d3heatmap(OverlapMatrix2, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }#first else if ends
      else if(col_include2 == length(colnames(OverlapTable3)) & gene_intersect2 == TRUE){
        OverlapMatrix2 <- data.matrix(OverlapTable3[!(gene_pos2), ])
        output$GeneHeatmap3 <- renderD3heatmap({
          d3heatmap(OverlapMatrix2, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }# second else if ends
      else if(col_include2 == length(colnames(OverlapTable3)) & gene_intersect2 == FALSE){
        OverlapMatrix2 <- data.matrix(OverlapTable3)
        output$GeneHeatmap3 <- renderD3heatmap({
          d3heatmap(OverlapMatrix2, Rowv = NA, Colv = NA, colors = brewer.pal(5, "Reds"),
                    scale = "row", xaxis_font_size=10, yaxis_font_size=10, na.rm = FALSE)
        })
      }
    }
  )
}
)
