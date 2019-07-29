voteTransformation = function(vote.0, vote.1){
final.vote = lapply(names(vote.0), function(i){
  # i = "2"
  cells.0 = as.character(vote.0[[i]]$cell)
  cells.1 = vote.1[as.integer(cells.0)]
  joint.data.frame = data.frame()
  joint.data.frame.1 = data.frame()
  for(j in 1:length(cells.1)){
    # j = 2
    temp.data.frame = data.frame("cell" = names(cells.1[[j]]), "simple.proportion" = vote.0[[i]]$simple.proportion[j]*as.vector(cells.1[[j]]),
                                 "compound.proportion" = vote.0[[i]]$compound.proportion[j]*as.vector(cells.1[[j]]))
    joint.data.frame = rbind(joint.data.frame, temp.data.frame)
  }
  cell.vote = joint.data.frame %>% group_by(cell) %>% summarise(simple.proportion = sum(simple.proportion), compound.proportion =
                                                                  sum(compound.proportion))
  cell.vote$simple.proportion = cell.vote$simple.proportion/sum(cell.vote$simple.proportion)
  cell.vote$compound.proportion = cell.vote$compound.proportion/sum(cell.vote$compound.proportion)
  cell.vote = as.data.frame(dplyr::arrange(cell.vote, desc(compound.proportion)))
  # print(cell.vote)
})
names(final.vote) = names(vote.0)
return(final.vote)
}
