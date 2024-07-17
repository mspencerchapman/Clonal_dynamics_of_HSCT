#----------------------------------
# Load packages (and install if they are not installed yet)
#----------------------------------
cran_packages=c("ggplot2","dplyr","RColorBrewer","tibble","stringr","readr","ggtext")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------
root_dir="~/R_work/Clonal_dynamics_of_HSCT"

#----------------------------------
# Set the ggplot2 theme for plotting
#----------------------------------

my_theme<-theme(text = element_text(family="Helvetica"),
                axis.text = element_text(size = 5),
                axis.title = element_text(size=7),
                legend.text = element_text(size=5),
                legend.title = element_text(size=7),
                strip.text = element_text(size=7),
                legend.spacing = unit(1,"mm"),
                legend.key.size= unit(5,"mm"))

#----------------------------------
# Read in the necessary metadata
#----------------------------------

Pair_metadata<-readr::read_csv(paste0(root_dir,"/data/metadata_files/Pair_metadata.csv"))
Pair_metadata$Pair_new<-factor(Pair_metadata$Pair_new,levels=paste("Pair",1:nrow(Pair_metadata),sep = "_"))
Pair_cols<-RColorBrewer::brewer.pal(10,"Paired")
names(Pair_cols)<-levels(Pair_metadata$Pair_new)

#----------------------------------
# Summary visualization of the time course
#----------------------------------

plot.timecourse<-Pair_metadata%>%
  ggplot(aes(x=Pair_new,xend=Pair_new,y = Age_at_transplant,yend=Age))+
  #geom_segment(size=1,linetype=1,col="darkred")+
  geom_segment(arrow=arrow(angle = 10,ends = "last",length=unit(0.4,"cm")),size=1,linetype=1,col="darkred")+
  geom_segment(aes(x=Pair_new,xend=Pair_new,y=0,yend=Age_at_transplant),col="darkgrey",lineend = "round")+
  geom_point(aes(x=Pair_new,y=0),col="darkgrey")+
  geom_point(aes(x=Pair_new,y=Age_at_transplant,shape=stem_cell_source),size=3,col="darkred")+
  scale_y_reverse(breaks = seq(0,80,10))+
  labs(x="",y="Age (years)", shape="Stem Cell\nSource")+
  theme_classic()+
  my_theme+
  theme(axis.line.x=element_blank(),axis.ticks.x = element_blank(),panel.grid.major.y = element_line())

#Visualization of the MNC & CD34+ dose where available
plot.cellcounts<-Pair_metadata%>%
  dplyr::select(Pair_new,MNC_dose,CD34_dose)%>%
  gather(-Pair_new,key="Cell_type",value="Cell_dose")%>%
  ggplot(aes(x=Pair_new,y=Cell_dose,fill=Cell_type))+
  geom_bar(stat="identity",position="dodge",col="black",size=0.3,width=0.5)+
  theme_bw()+
  scale_y_continuous(breaks=seq(0,20,2))+
  scale_fill_brewer(palette="Set2")+
  geom_hline(yintercept = 2,col=brewer.pal(3,"Set2")[1],linetype=2)+
  labs(x="Pair",y="Cell dose (CD34, x 10^(6)/kg; MNC, x 10^(8)/kg)",fill="Cell type")+
  my_theme+
  theme(axis.title.y=ggtext::element_markdown())
