# Wright-Fisher simulator written by Graham Coop

par(mar=c(2,2,2,2))

track_lineages<-function(N.vec, n.iter, num.tracked, col.allele=c("red","blue","orange","green","yellow","brown")){
num.gens<-length(N.vec)
offset<-0.2
for(iter in 1:n.iter){
	N.max<-max(N.vec)
	N<-N.vec[num.gens]
	N.prev<-N.vec[num.gens-1]
	plot(c(1,num.gens),c(1,N.max),type="n",axes=FALSE,xlab="",ylab="")
	mtext(side=1,line=1,"Generations")
	
	track.this.allele<-vector("list", 2*N)
	
	track.this.allele[sample(1:(2*N),num.tracked)]<-1:num.tracked
	
	track.this.allele.next.gen<-vector("list", 2*N.prev)
	
	for(i in num.gens:2){
		N<-N.vec[i]
		N.prev<-N.vec[i-1]
		track.this.allele.next.gen<-vector("list", 2*N.prev)
	for(ind in 1:N){
		
			par<-sample(1:N.prev,2,replace=FALSE)
			which.allele<-sample(c(-1,1),1)
			lines(c(i,i-1), c(ind-offset,par[1]+which.allele*offset),col="light grey",lwd=0.5)
			if(!is.null(track.this.allele[[2*ind-1]])){
				this.one<-2*par[1] +ifelse(which.allele==1,0,-1); 
				track.this.allele.next.gen[[this.one]]  <- c(track.this.allele.next.gen[[this.one]],track.this.allele[[2*ind-1]])
				}
				
			which.allele<-sample(c(-1,1),1)
			lines(c(i,i-1), c(ind+offset,par[2]+which.allele*offset),col="light grey",lwd=0.5)
			if(!is.null(track.this.allele[[2*ind]])){ 
						this.one<-2*par[2] +ifelse(which.allele==1,0,-1); 
				track.this.allele.next.gen[[ this.one]]  <- c(track.this.allele.next.gen[[this.one]],track.this.allele[[2*ind]])
				}
	#		recover()
		}
		for(this.allele in 1:num.tracked){ 
			daughter<-which(sapply(track.this.allele,function(allele){any(allele==this.allele)}))
			parent<-which(sapply(track.this.allele.next.gen,function(allele){any(allele==this.allele)}))
			lines(c(i,i-1), c(ceiling(daughter/2)+offset* ifelse(daughter %% 2,-1,1) ,ceiling(parent/2) + offset*ifelse(parent %% 2,-1,1) ),col=col.allele[this.allele],lwd=2)
			}
			
		points(rep(i,N),1:N+offset, pch=19,cex=1)
		points(rep(i,N),1:N-offset, pch=19,cex=1)
		 track.this.allele<-track.this.allele.next.gen
		}
		
	
 	}
	
}

