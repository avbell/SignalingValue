"bary.text" <-
function( point1, labels, lcex ) {
  ## default point is an empty circle
  pt <- bary.toscreen(point1[1], point1[2]);
  text( x=pt[1], y=pt[2], labels=labels, xpd=TRUE, cex=lcex )
}
