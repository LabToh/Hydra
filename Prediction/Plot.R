DiseasePlot = function(top1,top5,disease_name = "Meier-Gorlin syndrome",plot_disease_name = NULL,reported_gene_name=NULL){
	if(is.null(plot_disease_name)){
		plot_disease_name=disease_name
	}
	library(igraph)
	focus_disease_data_top1 = as.matrix(top1[top1[,2]==disease_name,])
	if(ncol(focus_disease_data_top1)==1){
		focus_disease_data_top1 = t(focus_disease_data_top1)
	}
	focus_disease_data = top5[top5[,2]==disease_name,]
	if(ncol(focus_disease_data)==1){
		focus_disease_data = t(focus_disease_data)
	}
	focus_disease_data[,2] = plot_disease_name
	net = graph(t(as.matrix(focus_disease_data)),directed=F)
	name_vector = V(net)
	pos = which(names(name_vector)==plot_disease_name)
	len = nrow(focus_disease_data)+1
	vertex_color = rep("white",len)
	vertex_color[pos] = "orange"
	vertex_label_color = rep("black",len)
	vertex_label_color[pos] = "white"
	vertex_label_cex = rep(0.5,len)
	vertex_label_cex[pos] = 0.7
	index = which(names(name_vector)%in%focus_disease_data_top1[,1])
	vertex_frame_color = rep("black",len)
	vertex_frame_color[index] = "red"
	vertex_shape = rep("circle",len)
	vertex_shape[pos] = "square"
	vertex_size = rep(25,len)
	vertex_size[pos] = 35
	if(!is.null(reported_gene_name)){
		site = which(names(name_vector)%in%reported_gene_name)
		vertex_color[site] = "blue"
		vertex_label_color[site] = "white"
	}
	cat("Top1 Genes = ",nrow(focus_disease_data_top1),": Top5 Genes = ",nrow(focus_disease_data),"\n")
	plot(net,edge.width=1,vertex.size=vertex_size,vertex.shape=vertex_shape,vertex.label.cex=vertex_label_cex,vertex.color=vertex_color,vertex.label.color=vertex_label_color,vertex.label.family="sans",vertex.frame.color=vertex_frame_color)
}

top1 = read.delim("/path/to/Prediction/eqhydraTop1PredictionPathwayIntegrated.txt",header=F)
top5 = read.delim("/path/to/Prediction/eqhydraTop5PredictionPathwayIntegrated.txt",header=F)

set.seed(127)
DiseasePlot(top1,top5,"Meier-Gorlin syndrome"," Meier-Gorlin \n syndrome",c("CDC45","MCM5"))
DiseasePlot(top1,top5,"Alzheimer's disease"," Alzheimer's \n disease")
