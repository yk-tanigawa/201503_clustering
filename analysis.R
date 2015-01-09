path = "~/moris/2_clustering/ytanigawa/"
fig_output_path = paste(path, "report/", sep="")

setwd(dir = path)
data <- read.table("results.txt")
name <-c("n", "d", "k", "time", "err", "rep")
names(data)=name
v_range<-rbind(c(1000,10000,100000), # 変数nの値域
               c(2,10,100),          # 変数dの値域
               c(10,100,1000))       # 変数kの値域

for(variable in 1:3){                # グラフの横軸にする変数
  for(observed in 4:6){              # グラフの縦軸にする変数
    v2<-((variable %% 3) + 1)        # グラフの点の形を変化させる変数
    v3<-(((variable + 1) %% 3) + 1)  # グラフの色を変化させる変数
    v2_range<-v_range[v2,]           # 変数の値域をセット
    v3_range<-v_range[v3,]           # 変数の値域をセット    
    png(paste(fig_output_path, name[variable],"_",name[observed],".png", sep = ""),
        pointsize = 30, width = 1200, height = 1200)
                                     # 保存するpng画像の設定    
    par(new=F)                       # 前の画像に重ね書きしない
    # par(xpd=T)                       # 凡例の枠外への描画を許可
    x_lim<-c(min(data[variable]), max(data[variable]))
    y_lim<-c(min(data[observed]), max(data[observed]))
                                     # グラフの描画範囲の設定
    plot(1,0, type="n", xlim=x_lim, ylim=y_lim, log="x", 
         main = paste(name[variable]," vs. ",name[observed], sep = ""),
         xlab = name[variable], ylab = name[observed])
                                     # 軸とタイトルのみからなるグラフを描画
    cols<-rainbow(3)
    for(v2_loop in 1:3){
      for(v3_loop in 1:3){
        par(new=T)
        plot_data<-data[intersect(which(data[v2]==v2_range[v2_loop]),
                                  which(data[v3]==v3_range[v3_loop])),]
                                     # 描画するデータを抽出
        plot(plot_data[,c(variable, observed)],
             xlim=x_lim, ylim=y_lim, ann=F, axes=F, log="x",
             pch=v2_loop, col=cols[v3_loop])
        lines(plot_data[,c(variable, observed)],
              xlim=x_lim, ylim=y_lim, ann=F,
              col=cols[v3_loop])
      }
    }
    legend("topleft",bty = "n",inset = c(0, 0),
           legend = c(paste(name[v2], " = ", v2_range[1],sep=""),
                      paste(name[v2], " = ", v2_range[2],sep=""),
                      paste(name[v2], " = ", v2_range[3],sep="")),
           pch = 1:3)
    legend("topleft",bty = "n",inset = c(0, 0.15),
           legend = c(paste(name[v3], " = ", v3_range[1],sep=""),
                      paste(name[v3], " = ", v3_range[2],sep=""),
                      paste(name[v3], " = ", v3_range[3],sep="")),
           lty = 1, col = cols)
    dev.off()
  }
}
