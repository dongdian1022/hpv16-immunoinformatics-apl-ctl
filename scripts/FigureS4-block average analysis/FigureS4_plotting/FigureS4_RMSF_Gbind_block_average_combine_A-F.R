# ==========================================
# Figure S4
#
# Block-averaged MD robustness analysis
#
# A-C: RMSF block averaging
# D-F: MM-GBSA block averaging
#
# Final v7
#
# ==========================================


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)



# ==========================================
# Working directory
# ==========================================

setwd("E:/HPV-T非临床文章/0.文章投稿+确定版本/vaccines投稿后返修等工作/gpt返修/1.block average/Figure7_plot/summary")



# ==========================================
# Theme
# ==========================================

theme_S4 <- theme_classic(
  base_size = 14
)+
  theme(
    
    axis.line =
      element_line(
        color="black",
        linewidth=1
      ),
    
    axis.text.x =
      element_text(
        color="black",
        size=11,
        face="bold"
      ),
    
    axis.text.y =
      element_text(
        color="black",
        size=11,
        face="bold"
      ),
    
    axis.title.x =
      element_text(
        color="black",
        size=13,
        face="bold"
      ),
    
    axis.title.y =
      element_text(
        color="black",
        size=13,
        face="bold"
      ),
    
    axis.ticks =
      element_line(
        color="black",
        linewidth=1
      ),
    
    axis.ticks.length =
      unit(
        0.18,
        "cm"
      ),
    
    legend.position =
      "top",
    
    legend.title =
      element_blank(),
    
    legend.text =
      element_text(
        size=10,
        face="bold"
      ),
    
    legend.key.width =
      unit(
        1.0,
        "cm"
      ),
    
    plot.margin =
      margin(
        8,8,8,8,
        unit="pt"
      )
    
  )



# ==========================================
# Read RMSF
# ==========================================

read_rmsf <- function(
    file,
    group
){
  
  dat <- read.table(
    file,
    header=FALSE,
    comment.char="#"
  )
  
  
  colnames(dat)<-
    c(
      "Position",
      "RMSF",
      "SEM"
    )
  
  
  dat$Group <- group
  
  
  return(dat)
  
}




# ==========================================
# Mutation annotation
# ==========================================

add_mutation <- function(
    p,
    x,
    label
){
  
  
  for(i in seq_along(x)){
    
    
    p <- p +
      
      annotate(
        "text",
        x=x[i],
        y=2.25,
        label=label[i],
        size=3.8,
        fontface="bold"
      )+
      
      
      annotate(
        "segment",
        x=x[i],
        xend=x[i],
        y=2.12,
        yend=1.95,
        linewidth=0.7,
        arrow=
          arrow(
            length=
              unit(
                0.1,
                "inches"
              ),
            type="closed"
          )
      )
    
  }
  
  
  return(p)
  
}





# ==========================================
# RMSF plot
# ==========================================

make_rmsf <- function(
    
  WT_file,
  APL_file,
  
  wt_label,
  apl_label,
  
  wt_color,
  apl_color,
  
  peptide_length,
  
  anchor,
  
  mutation_x,
  mutation_label
  
){
  
  
  wt <- read_rmsf(
    WT_file,
    "WT"
  )
  
  
  apl <- read_rmsf(
    APL_file,
    "APL"
  )
  
  
  
  df <- bind_rows(
    wt,
    apl
  )
  
  
  
  df$Group <- factor(
    df$Group,
    levels=c(
      "APL",
      "WT"
    )
  )
  
  
  
  
  p <- ggplot(
    
    df,
    
    aes(
      x=Position,
      y=RMSF,
      color=Group,
      group=Group
    )
    
  )+
    
    
    annotate(
      "rect",
      xmin=anchor[1]-0.45,
      xmax=anchor[1]+0.45,
      ymin=0,
      ymax=2.5,
      fill="grey80",
      alpha=0.6
    )+
    
    
    annotate(
      "rect",
      xmin=anchor[2]-0.45,
      xmax=anchor[2]+0.45,
      ymin=0,
      ymax=2.5,
      fill="grey80",
      alpha=0.6
    )+
    
    
    geom_errorbar(
      aes(
        ymin=RMSF-SEM,
        ymax=RMSF+SEM
      ),
      width=0.1,
      linewidth=0.5
    )+
    
    
    geom_line(
      linewidth=1.1
    )+
    
    
    geom_point(
      size=2.8
    )+
    
    
    scale_color_manual(
      
      values=c(
        APL=apl_color,
        WT=wt_color
      ),
      
      labels=c(
        APL=apl_label,
        WT=wt_label
      )
      
    )+
    
    
    scale_x_continuous(
      
      breaks=1:peptide_length,
      
      labels=paste0(
        "P",
        1:peptide_length
      ),
      
      limits=c(
        0.5,
        peptide_length+0.5
      )
      
    )+
    
    
    scale_y_continuous(
      
      limits=c(
        0,
        2.5
      ),
      
      breaks=seq(
        0,
        2.5,
        0.5
      )
      
    )+
    
    
    labs(
      
      x="Peptide position",
      
      y="Backbone RMSF (Å)"
      
    )+
    
    
    theme_S4+
    
    
    guides(
      
      color=
        guide_legend(
          nrow=1,
          override.aes=list(
            size=3,
            linewidth=1.2
          )
        )
      
    )
  
  
  
  p <- add_mutation(
    p,
    mutation_x,
    mutation_label
  )
  
  
  
  return(p)
  
}





# ==========================================
# MMGBSA plot
# ==========================================

make_mmgbsa <- function(
    
  WT,
  APL,
  
  wt_label,
  apl_label,
  
  wt_color,
  apl_color,
  
  ylim
  
){
  
  
  
  df <- data.frame(
    
    Block=factor(
      
      c(
        "20–30 ns",
        "30–40 ns",
        "40–50 ns"
      ),
      
      levels=c(
        "20–30 ns",
        "30–40 ns",
        "40–50 ns"
      )
      
    ),
    
    WT=WT,
    
    APL=APL
    
  )
  
  
  
  df <- df %>%
    
    pivot_longer(
      
      cols=c(
        WT,
        APL
      ),
      
      names_to="Group",
      
      values_to="dG"
      
    )
  
  
  
  df$Group <- factor(
    
    df$Group,
    
    levels=c(
      "APL",
      "WT"
    )
    
  )
  
  
  
  
  ggplot(
    
    df,
    
    aes(
      x=Block,
      y=dG,
      color=Group,
      group=Group
    )
    
  )+
    
    
    geom_line(
      linewidth=1.2
    )+
    
    
    geom_point(
      size=3
    )+
    
    
    scale_color_manual(
      
      values=c(
        APL=apl_color,
        WT=wt_color
      ),
      
      labels=c(
        APL=apl_label,
        WT=wt_label
      )
      
    )+
    
    
    scale_y_continuous(
      
      limits=ylim,
      
      breaks=seq(
        ylim[1],
        ylim[2],
        5
      )
      
    )+
    
    
    labs(
      
      x="Trajectory block",
      
      y="ΔGbind (kcal/mol)"
      
    )+
    
    
    theme_S4+
    
    
    guides(
      
      color=
        guide_legend(
          
          nrow=1,
          
          override.aes=list(
            size=3,
            linewidth=1.2
          )
          
        )
      
    )
  
  
}





# ==========================================
# A-C RMSF
# ==========================================


A <- make_rmsf(
  
  "E6wt_10mer_rmsf_meanSEM.dat",
  
  "E6apl_10mer_rmsf_meanSEM.dat",
  
  "WT (E6wt)",
  
  "APL (E6apl)",
  
  "#1A4F8A",
  
  "#8B3A0F",
  
  10,
  
  c(2,10),
  
  c(1,2),
  
  c("T1Y","I2L")
  
)



B <- make_rmsf(
  
  "E7wt1_9mer_rmsf_meanSEM.dat",
  
  "E7apl1_9mer_rmsf_meanSEM.dat",
  
  "WT (E7wt1)",
  
  "APL (E7apl1)",
  
  "#5BA3D0",
  
  "#E08A3C",
  
  9,
  
  c(2,9),
  
  c(1),
  
  c("L1Y")
  
)



C <- make_rmsf(
  
  "E7wt2_10mer_rmsf_meanSEM.dat",
  
  "E7apl2_10mer_rmsf_meanSEM.dat",
  
  "WT (E7wt2)",
  
  "APL (E7apl2)",
  
  "#3A7EBF",
  
  "#C1622A",
  
  10,
  
  c(2,10),
  
  c(2,10),
  
  c("M2L","T10V")
  
)





# ==========================================
# D-F MMGBSA
# ==========================================


D <- make_mmgbsa(
  
  c(-90.07,-85.64,-76.71),
  
  c(-69.21,-72.28,-79.28),
  
  "WT (E6wt)",
  
  "APL (E6apl)",
  
  "#1A4F8A",
  
  "#8B3A0F",
  
  c(-95,-65)
  
)



E <- make_mmgbsa(
  
  c(-76.91,-80.50,-68.98),
  
  c(-84.40,-76.79,-79.45),
  
  "WT (E7wt1)",
  
  "APL (E7apl1)",
  
  "#5BA3D0",
  
  "#E08A3C",
  
  c(-90,-65)
  
)



F <- make_mmgbsa(
  
  c(-83.39,-79.38,-71.51),
  
  c(-86.71,-84.54,-75.45),
  
  "WT (E7wt2)",
  
  "APL (E7apl2)",
  
  "#3A7EBF",
  
  "#C1622A",
  
  c(-90,-70)
  
)




# ==========================================
# Combine
# ==========================================


FigureS4 <-
  
  (A+B+C)/(D+E+F)+
  
  plot_annotation(
    
    tag_levels="A",
    
    theme = theme(
      
      plot.tag =
        element_text(
          size=16,
          face="bold",
          color="black"
        ),
      
      plot.tag.position =
        "topleft"
      
    )
    
  )



FigureS4



# ==========================================
# Export
# ==========================================


ggsave(
  
  "FigureS4_block_average_final_v7.pdf",
  
  FigureS4,
  
  width=12,
  
  height=8.5
  
)



ggsave(
  
  "FigureS4_block_average_final_v7.tiff",
  
  FigureS4,
  
  width=12,
  
  height=8.5,
  
  dpi=600,
  
  compression="lzw"
  
)



ggsave(
  
  "FigureS4_block_average_final_v7.png",
  
  FigureS4,
  
  width=12,
  
  height=8.5,
  
  dpi=600
  
)