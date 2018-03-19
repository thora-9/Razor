euler <- function(f, y0, a, b, h)
{
  t <- a
  y <- y0
  
  while (t < b)
  {
    cat(sprintf("%6.3f %6.3f\n", t, y))
    t <- t + h
    y <- y + h*f(t, y)
  }
}