#include "mylib.h"

string int2str(int i) {
        if (i<0) {
                return "-" + int2str(-i);
        } else if (i<10) {
                return i2s[i];
        } else {
                return int2str(i/10)+i2s[i%10];
        }
}

void tokenizer(string seq, string separators, vector<string>* result) {
    // split the seq into many parts by "separators"
    // the vector<string> *result cannot be NULL
    result->clear();
    int startpos = (int) seq.find_first_not_of(separators);
    while (startpos != (int) string::npos) {
        int endpos = (int) seq.find_first_of(separators, startpos);
        if (endpos != (int) string::npos) {
            result->push_back(seq.substr(startpos, endpos-startpos));
            startpos = (int) seq.find_first_not_of(separators, endpos);
        } else {
            result->push_back(seq.substr(startpos));
            break;
        }
    }
}

void randomizePos(vector<int>& arr, int n) {
        // reorder the items inside the array randomly and only report the top n items
        int i,r,t;
        for (i=0; i<n; i++) {
                r=rand() % (arr.size()-i);
                // swap between pos i and pos r+i
                t = arr[i];
                arr[i] = arr[r+i];
                arr[r+i] = t;
        }
}

void randomizePos(vector<int>& arr) {
        // reorder the items inside the array randomly
	randomizePos(arr, arr.size()-1);
}

void statistics(long double* arr, int n, long double& min, long double& max, long double& avg, long double& median, long double& stdev) {
	// get the minimum, maximum, average, median, and standard deviation for the set of numbers
	vector<long double> a;
	int i;
	for (i=0; i<n; i++)
		a.push_back(arr[i]);
	statistics(a, min, max, avg, median, stdev);
}

void statistics(vector<long double>& arr, long double& min, long double& max, long double& avg, long double& median, long double& stdev) {
	// get the minimum, maximum, average, median, and standard deviation for the set of numbers
	if (arr.size() == 0)
		return;
	int i,n;
	long double t;
	sort(arr.begin(), arr.end());
	long double sum = arr[0];
	min = arr[0];
	max = arr[0];
	n = arr.size();
	for (i=1; i<n; i++) {
		sum += arr[i];
		if (arr[i] < min)
			min = arr[i];
		else if (arr[i] > max)
			max = arr[i];
	}
	avg = sum / n;
	if (n%2 == 1)
		median = arr[(n-1)/2];
	else
		median = (arr[n/2]+arr[n/2-1])/2.0;
	stdev = 0.0;
	for (i=0; i<n; i++) {
		t = arr[i]-avg;
		stdev += (t*t);
	}
	stdev = sqrtl(stdev/(n-1));
}

int getTopK(vector<int>& arr, int topK) {
	// return the value of topK
	vector<int> t;
	int i;
	for (i=0; i<arr.size(); i++)
		t.push_back(arr[i]);
	sort(t.begin(), t.end());
	return t[topK];
}

long subcombin(int n, int k, int dim, long* ans) {
	if (ans[n*dim + k] == 0)
		ans[n*dim + k] = subcombin(n-1, k-1, dim, ans) + subcombin(n-1, k, dim, ans);
	return ans[n*dim + k];
}

long double logfact(int n) {
	
	if (initialCumLogI && n < MAXDIM)
		return cumlogI[n];
		
	long double i = 0.0;
	for (int k=2; k<=n; k++) {
		i += logl((long double)k);
	}
	return i;
}

long double logfact(int n, int k) {
	
	if (initialCumLogI==1 && n < MAXDIM)
		return cumlogI[n] - cumlogI[k-1];
		
	long double i = 0.0;
	for (int m=k; m<=n; m++) {
		i += logl((long double)m);
	}
	return i;
}

long double logCombin(int n, int k) {
	// = log(n!)-log(k!)-log((n-k)!)
	if (n==k)
		return 0.0;
	return logfact(n,k+1)-logfact(n-k);
}

long double logBinomial(int x, int n, long double p) {
	return logCombin(n,x) + x*logl(p) + (n-x)*logl(1.0-p);
	// return x*logl(p) + (n-x)*logl(1.0-p);
}


long combin(int n, int k) {
	// n >= k
	int i,j;
	if (k>n || n==0)
		return 0;
	if (k==1)
		return n;
	if (k==0 || k==n)
		return 1;
	long *ans = new long[n*k];
	memset(ans, 0, n*k*sizeof(long));
	for (i=0; i<=k; i++)
		ans[i*k+i] = 1;
	for (i=0; i<n; i++)
		ans[i*k] = i+1;
	ans[0] = 1;
	
	long out = subcombin(n-1, k-1, k, ans);
	delete[] ans;
	
	return out;
}

void removeSpaces(string& s) {
	// remove all spaces inside the string
	int i,k;
	k=0;
	for (i=0; i<s.length(); i++) {
		if (s[i]!=' ') {
			if (k<i)
				s[k]=s[i];
			k++;
		}
	}
	if (k==0)
		s = "";
	else if (k<s.length())
		s.resize(k);
}

// compute the BIC value
long double BIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl ((long double) dataSize);
}

// compute the AIC value
long double AIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * 2.0;
}

// compute the Adjusted BIC value
long double ABIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl (((long double) dataSize + 2.0)/24.0);
}


// integer to string
// minLen : minimum length of the integer.
// If the number of digits < minimum length, then the integer will be displayed with leading zeros, unless the integer is zero
string intToStr(int i, int minLen) {
	char* d2c = (char*) "0123456789"; // the array for digit to char
	if (i<0) {
		return "-" + intToStr(-i, minLen);
	} else if (i<10) {
		if (minLen > 0)
			return string(minLen-1, '0') + string(1,d2c[i]);
		else
			return string(1,d2c[i]);
	} else {
		return intToStr(i/10, minLen-1)+string(1,d2c[i%10]);
	}
}

// double to string
string doublToStr(double d, double decimalPlace) {
	if (isnan(d)) {
		return "N/A";
	} else if (d < 0) {
		return "-" + doublToStr(-d, decimalPlace);
	} else if (decimalPlace > 0) {
		d = roundNumber(d, decimalPlace);
		return int2str((int)d) + "." + intToStr((int)round((d - (int)d)*pow(10.0,decimalPlace)), decimalPlace);
	} else {
		return int2str((int)d);
	}
}

// round the double to a certain number of decimal place
double roundNumber(double d, int decimalPlace) {
	return round( d * pow(10.0, decimalPlace) ) / pow(10.0, decimalPlace);
}


// compute the CAIC value
long double CAIC(long double logL, int degFree, int dataSize) {
	return -2.0 * logL + (long double) degFree * logl ((long double) dataSize + 1.0);
}

// compute the AICc value
long double AICc(long double logL, int degFree, int dataSize) {
	return AIC(logL, degFree, dataSize) + (long double) 2.0 * degFree * (degFree + 1) / (dataSize - degFree - 1);
}

// initialize the cumlogI array
void initializeCumLogI() {
	int i;
	cumlogI[1] = 0;
	for (i=2; i<MAXDIM; i++) {
		cumlogI[i] = logl((long double)i) + cumlogI[i-1];
	}
	initialCumLogI = true;
}

// compute log(x+y), given log(x) and log(y)
long double log_x_plus_y(long double logx, long double logy) {
	long double loga, logb;
	if (logx <= logy) {
		loga = logx;
		logb = logy;
	} else {
		loga = logy;
		logb = logx;
	}
	return logl(expl(loga - logb) + 1) + logb;
}

