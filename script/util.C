double pull(double x, double mean, double err){

	double diff = x - mean;
	
	return diff * diff / ( 2 * err * err);

}

double min(double a, double b){
	if(a > b)
		return b;
	else
		return a;
}

double max(double a, double b){
	if(a > b)
		return a;
	else
		return b;
}

double sigmoid(double x){
	if(x > 0)
		return 1 / ( 1 + exp(-x));
	else{
		double _tmp = exp(x);
		return _tmp / ( 1 + _tmp );
	}
}

double neumaierSum( double *vec, size_t size ){
	double sum = vec[0];
	double c = 0;
	for(size_t i=1;i<size;++i){
		double t = sum + vec[i];
		if( fabs(sum) >= fabs(vec[i]) )
			c += (sum - t) + vec[i];
		else
			c += (vec[i] - t) + sum;
		sum = t;
	}
	return sum + c;
}


