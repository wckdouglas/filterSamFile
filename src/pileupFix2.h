#include <string>
#include <map>
#include <sstream>


using namespace std;
typedef map<string, int> baseCounter;

int countDigits(int number) 
{
	if (number < 10) 
	{
		return 1;
	}
	int count = 0;
	while (number > 0) 
	{
		number /= 10;
		count++;
	}
	return count;
}

string reverseComplement(string base){
    string out = base ;
    if (base == "A")
    {
        out = "t";
    }
    else if (base == "C")
    {
        out = "g";
    }
    else if (base == "G")
    {
        out = "c";
    }
    else if (base == "T")
    {
        out = "a";
    }
    return out;
}

void baseCount(baseCounter &counter, string ref, char readPos)
{
    string base;
    stringstream ss;
    ss << readPos;
    ss >> base; 
    if (base == ".")
    {
        base = ref;
    }
    else if (base == ",")
    {
        base = reverseComplement(ref);
    }
    counter[base] ++;
}

void fixpileup(baseCounter &counter,
				int &deletion, int &insertion, 
				string reads, string baseQuals, string ref, 
				int qualThreshold, int &cov,
				int &start, int &end)
{
	int i = 0, j = 0, current = 0, qual;
	char readPos;
    while (i < reads.length())
    {
		readPos = reads.at(i);
        if (readPos == '+')
		//insertion
        {
			i ++ ; 
			current = 0;
			insertion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			i += current - countDigits(current);
        }
		else if (readPos == '-')
		// deletion
		{
			i ++ ; 
			current = 0;
			deletion ++;
			while (isdigit(reads.at(i)))
			{
				current += current * 10 + (reads[i]-'0');
				i++;
			}
			i += current - countDigits(current);
		}
        else if (readPos == '^')
        {
            i ++;
			if (islower(reads.at(i+1)) || reads.at(i+1) == ',')
			{
				start ++;
			}
			else if (isupper(reads.at(i+1)) || reads.at(i+1) == '.')
			{
				end ++;
			}
        }
        else if (readPos == '$')
        {
			if (islower(reads.at(i-1)) || reads.at(i-1) == ',')
			{
				end ++;
			}
			else if (isupper(reads.at(i-1)) || reads.at(i-1) == '.')
			{
				start ++;
			}
        }
		else if (readPos == '<' || readPos == '>')
		{
			cov --;
		}
		else 
		{
            qual = baseQuals[j] - 33 ;
			j++;
			if (qual < qualThreshold || readPos == '*')
			{
				cov --;
			}
			else 
			{
				baseCount(counter, ref, readPos);
			}
		}
		i++;
    }
}
