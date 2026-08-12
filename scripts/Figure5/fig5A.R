# Figure5A_v8_final_Welch_ttest.R
# Based on original Figure5A_v5.R
# Modification: add WT-APL statistical comparisons only.
# Original plotting style preserved:
# - no floating bars
# - original WT/APL legend
# - original colors and axis settings

library(ggplot2)
library(dplyr)
library(grid)

# ===============================
# Original data and plotting code
# ===============================

pal <- c(
  "E7wt2"="#3A7EBF",
  "E7apl2"="#C1622A",
  "E7apl3"="#D4732A",
  "E6wt"="#1A4F8A",
  "E6apl"="#8B3A0F",
  "Pos Ctrl"="#7A7A7A"
)

col_wt <- "#3A7EBF"
col_apl <- "#C1622A"

df_raw <- data.frame(
  peptide=rep(c("E7wt2","E7apl2","E7apl3","E6wt","E6apl","Pos Ctrl"),each=2),
  index=c(1.392,1.545,3.181,2.901,2.239,2.674,
          1.477,1.348,2.555,2.341,2.840,3.130)
)

pep_order <- c("E7wt2","E7apl2","E7apl3","E6wt","E6apl","Pos Ctrl")

df_b <- data.frame(
  peptide=pep_order,
  group=factor(c("WT","APL","APL","WT","APL","CTRL"),
               levels=c("WT","APL","CTRL")),
  xpos=c(1.0,1.7,2.4,4.0,4.7,6.3),
  mean=sapply(pep_order,function(x) mean(df_raw$index[df_raw$peptide==x])),
  sd=sapply(pep_order,function(x) sd(df_raw$index[df_raw$peptide==x])),
  bar_col=unname(pal[pep_order])
)

legend_df <- data.frame(
  group=factor(c("WT","APL"),levels=c("WT","APL")),
  x=c(-99,-99),
  y=c(-99,-99)
)

# ===============================
# Welch t tests
# ===============================

p_E7apl2 <- t.test(index~peptide,
                   data=subset(df_raw, peptide %in% c("E7wt2","E7apl2")),
                   var.equal=FALSE)$p.value

p_E7apl3 <- t.test(index~peptide,
                   data=subset(df_raw, peptide %in% c("E7wt2","E7apl3")),
                   var.equal=FALSE)$p.value

p_E6 <- t.test(index~peptide,
               data=subset(df_raw, peptide %in% c("E6wt","E6apl")),
               var.equal=FALSE)$p.value

# ===============================
# Significance brackets
# ===============================

bracket <- data.frame(
  x1=c(1.0,1.0,4.0),
  x2=c(1.7,2.4,4.7),
  y=c(3.55,3.82,3.35),
  label=c("*","*","*")
)

# ===============================
# Plot
# ===============================

p <- ggplot()+
  geom_hline(yintercept=1,linetype="dashed",
             color="#9A9A9A",linewidth=0.5)+
  geom_bar(data=df_b,aes(x=xpos,y=mean),
           stat="identity",width=0.60,
           fill=df_b$bar_col,color=NA)+
  geom_errorbar(data=df_b,
                aes(x=xpos,ymin=mean-sd,ymax=mean+sd),
                width=0.18,linewidth=0.6)+
  geom_text(data=df_b,
            aes(x=xpos,y=mean+sd+0.12,label=sprintf("%.1f",mean)),
            size=3,fontface="bold")+
  geom_bar(data=legend_df,
           aes(x=x,y=y,fill=group),
           stat="identity")+
  geom_segment(data=bracket,
               aes(x=x1,xend=x2,y=y,yend=y),
               linewidth=0.5)+
  geom_segment(data=bracket,
               aes(x=x1,xend=x1,y=y-0.08,yend=y),
               linewidth=0.5)+
  geom_segment(data=bracket,
               aes(x=x2,xend=x2,y=y-0.08,yend=y),
               linewidth=0.5)+
  geom_text(data=bracket,
            aes(x=(x1+x2)/2,y=y+0.12,label=label),
            size=4)+
  scale_fill_manual(values=c(WT=col_wt,APL=col_apl),
                    name=NULL)+
  scale_x_continuous(
    breaks=df_b$xpos,
    labels=df_b$peptide,
    expand=expansion(add=c(0.5,0.5))
  )+
  scale_y_continuous(
    limits=c(0,3.8),
    breaks=c(0,1,2,3),
    expand=c(0,0)
  )+
  labs(x=NULL,y="MFI enhancement index")+
  theme_classic(base_size=12)+
  theme(
    axis.title.y=element_text(size=12,face="bold"),
    legend.position=c(0.88,1.01),
    legend.justification=c(1,1),
    legend.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    legend.key = element_rect(
      fill = "transparent",
      color = NA
    )
  )

print(p)

ggsave("Figure5A_v8_final.png",p,width=5.8,height=4.4,dpi=600,bg="white")
ggsave("Figure5A_v8_final.pdf",p,width=5.8,height=4.4,bg="white")
ggsave("Figure5A_v8_final.tiff",p,width=5.8,height=4.4,dpi=600,
       compression="lzw",bg="white")

message("Done!")