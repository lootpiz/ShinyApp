library(ggplot2)
library(shiny)
library(NMF)

data <- readRDS(file.path("./data", "TCGA-GBM_Expression_DESEQ2.RDS"))
gtex <- readRDS(file.path("./data", "GTEx-Cortex_Expression_DESEQ2.RDS"))
nmfr <- readRDS(file.path("./data", "TCGA_NMFresults.RDS"))

rank <- 5

ui <- fluidPage(
        titlePanel("TCGA-GBM Bulk RNA-seq"),
        sidebarLayout(
                sidebarPanel(
                        helpText("A combined cohort of TCGA-GBM and GTEx samples downloaded from Xenabrowser. Gene expression profile is RSEM expected_count output normalized using DESeq2."),
                        textInput(
                                inputId = "gene",
                                label = h3("Gene symbol"),
                                value = "TNC"
                        ),
                        downloadLink(
                                outputId = "export",
                                label = "Export Figure(s)"
                        )       
                ),
                mainPanel(
                        h2("TCGA-GBM, N=164 vs. GTEx-Cortex, N=105"),
                        tabsetPanel(
                                type = "tabs",
                                tabPanel("NMF subgroups", h5("Subgroups were defined via NMF (Non-negatvie Matrix Factorization)"), plotOutput("violin")),
                                tabPanel("GBM vs. Brain", plotOutput("violin2"))
                        )                        
                )
        )
)

server <- function(input, output, session) {
        generateViolinNMF <- function() {
                gene <- input$gene
                if (length(which(rownames(gtex) == gene)) > 0) {
                        results <- data.frame(
                                Group = "GTEx",
                                Expression = as.numeric(gtex[which(rownames(gtex) == gene), ])
                        )

                        subgroups <- predict(get(as.character(rank), nmfr$fit), "consensus")

                        for (subGroupLv in c(1:rank)) {
                                dataSubIdx <- which(colnames(data) %in% names(subgroups)[which(subgroups == subGroupLv)])

                                df <- data.frame(
                                        Group = paste0("Group_", subGroupLv),
                                        Expression = as.numeric(unlist(data[which(rownames(data) == gene), dataSubIdx]))
                                )

                                results <- rbind(results, df)
                        }

                        results$Group <- factor(results$Group, levels = c(paste0("Group_", c(1:rank)), "GTEx"))

                        dodge <- position_dodge(width = 0.5)
                        ggplot(data = results, aes(x = Group, y = Expression, fill = Group)) +
                                geom_violin(position = dodge, size = 0) +
                                geom_boxplot(width = 0.1, position = dodge, fill = "white") +
                                scale_fill_manual(values = c("gold", "orange", "brown1", "brown4", "indianred3", "grey80")) +
                                labs(
                                        title = gene,
                                        x = "", y = "Expression, Log2"
                                ) +
                                theme_bw() +
                                theme(
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                                        legend.position = "none",
                                        text = element_text(size = 20)
                                )
                }
        }
        generateViolinTvsB <- function() {
                gene <- input$gene
                if (length(which(rownames(gtex) == gene)) > 0) {
                        resultsGtex <- data.frame(
                                Group = "GTEx-Cortex",
                                Expression = as.numeric(gtex[which(rownames(gtex) == gene), ])
                        )

                        resultsGbm <- data.frame(
                                Group = "TCGA-GBM",
                                Expression = as.numeric(data[which(rownames(data) == gene), ])
                        )

                        results <- rbind(resultsGbm, resultsGtex)
                        results$Group <- factor(results$Group, levels = c("TCGA-GBM", "GTEx-Cortex"))

                        dodge <- position_dodge(width = 0.5)
                        ggplot(data = results, aes(x = Group, y = Expression, fill = Group)) +
                                geom_violin(position = dodge, size = 0) +
                                geom_boxplot(width = 0.1, position = dodge, fill = "white") +
                                scale_fill_manual(values = c("chartreuse2", "grey80")) +
                                labs(
                                        title = gene,
                                        x = "", y = "Expression, Log2"
                                ) +
                                theme_bw() +
                                theme(
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                                        legend.position = "none",
                                        text = element_text(size = 20)
                                )
                }
        }
        output$violin <- renderPlot({ generateViolinNMF() }, height = 600, width = 800)

        output$violin2 <- renderPlot({ generateViolinTvsB() }, height = 600, width = 800)

        output$export <- downloadHandler(
                filename = function() {
                        paste0("TCGA_", input$gene, "_", Sys.Date(), ".pdf")
                },
                content = function(file) {
                        pdf(file = file)
                        print(generateViolinNMF())
                        print(generateViolinTvsB())
                        dev.off()
                }
        )
}

shinyApp(ui = ui, server = server)
