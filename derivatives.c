

double fp ( double ip, double im, double h )
{
  return (ip - im)/(2*h);
}

double fpp (double im, double i, double ip, double h)
{
  double der1 = fp(ip,i ,h/2);
  double der2 = fp(i ,im,h/2);
  return (der1  - der2)/h;
}

double fppp(double ip, double im, double i, double ipp, double imm, double h)
{
  double derr1 = fpp(ipp,i,ip,h);
  double derr2 = fpp(i,imm,im,h);
  return (derr1 - derr2)/(2*h);
}

