#define _USE_MATH_DEFINES

#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include "mpi.h"
#include <ctime>
#include <omp.h>

using namespace std;

//this is a function to generate a linspace of points
vector<double> linspace(double start, double end, int numPoints);
//function generates a vector of equally spaced points
vector<double> equalSpacing(double start, double spacing, int numPoints);
//generates random gaussian numbres
double gaus_box_mul_gen();
//price for the monte carlo call
double monte_call(int number_sims, double underlying_price, double strike_price, double risk_free, double volatility, double time);
//price for the monte carlo put
double monte_put(int num_sims, double underlying_price, double strike_price, double risk_free, double volatility, double time);
//formula for pulling from normal cdf
double normal_cdf(double x);
//this is for calcualting d_1 and d_2 for the sto proceess
double d_j(int which_one, double underlying_price, double strike_price, double risk_free, double volatility, double time);
//sto process call price validation
double call_price(double underlying_price, double strike_price, double risk_free, double volatility, double time);
//sto process put price validation
double put_price(double underlying_price, double strike_price, double risk_free, double volatility, double time);


//this function is used to generate a linspace over a given span
vector<double> linspace(double start, double end, int numPoints)
{
    //a vector is declared to store the values
    vector<double> ret_array;
    //after this the step size is calculated
    double step = (end - start) / (numPoints - 1);
    //after this a while loop is used while the point your currently at is less than the end
    while (start <= end) 
    {
        //the value is pushed onto the back of the vecotr
        ret_array.push_back(start);
        //the current position is incremented by the step size
        start += step;
    }
    //the vector is returned
    return ret_array;
}
//this creates a vector of equal spacing between each point
vector<double> equalSpacing(double start, double spacing, int numPoints) 
{
    //the return array is declared
    vector<double> ret_array;
    //using a for loop based off the number of points values are pushed onto the vector using the start plus the index times spacing
    for (int i = 0; i < numPoints; i++) 
    {
        //values are pushed on here
        ret_array.push_back(start + spacing*i);
    }
    //vector is returned
    return ret_array;
}

//implementation of box-muller algorithm, it is used to generate gaussian rnadom numbers
double gaus_box_mul_gen() 
{
    //declaring the basis for method
    //x declaration
    double x = 0.0;
    //y declaration
    double y = 0.0;
    //declaring euclidan square
    double euclid_sq = 0.0;

    
    //using do loop to generate two uniform random variables until the square of their distance is less than 1
    do {
        //assigning x value based off rands and max rand
        x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        //assigning y to be the same
        y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        //assigning the value of the square based off of x and y
        euclid_sq = x*x + y*y;
        //checking the condition for exit
    } while (euclid_sq >= 1.0);
    //returns x*(-2*log(euclid_sq)/euclid_sq)^.5
    return x*sqrt(-2 * log(euclid_sq) / euclid_sq);
}


//pricing method for each individual call for the monte carlo method
double monte_call(int number_sims, double underlying_price, double strike_price, double risk_free, double volatility, double time) 
{
    //declaring a value for the price adjusted based off the underling and the formula outlined in paper
    double price_adjusted = underlying_price * exp(time*(risk_free - 0.5*volatility*volatility));
    //setting current equal to 0 as starting point
    double price_current = 0.0;
    //initalizing the payout sum
    double payoff_total = 0.0;
    //for loop to go through all simultions
    for (int i = 0; i<number_sims; i++)
    {
        //declare a holder for the values form the random value generator
        double gauss_bm = gaus_box_mul_gen();
        //assign the current price based off calculation
        price_current = price_adjusted * exp(sqrt(volatility*volatility*time)*gauss_bm);
        //increasing the payoff total by what the max of the price is at that time
        payoff_total += max(price_current - strike_price, 0.0);
    }
    //return the average of the payoffs in effect coupling with the risk free return metric
    return (payoff_total / static_cast<double>(number_sims)) * exp(-risk_free*time);
}

//pricing method for each individual put for the monte carlo method
double monte_put(int num_sims, double underlying_price, double strike_price, double risk_free, double volatility, double time)
{
    //declaring a value for the price adjusted based off the underling and the formula outlined in paper\
    //using the put call parity when comparing put versus call
    double S_adjust = underlying_price * exp(time*(risk_free - 0.5*volatility*volatility));
    //setting current equal to 0 as starting point
    double S_cur = 0.0;
    //initalizing the payout sum
    double payoff_sum = 0.0;
    //for loop to go through all simultions
    for (int i = 0; i<num_sims; i++)
    {
        //declare a holder for the values form the random value generator
        double gauss_bm = gaus_box_mul_gen();
        //assign the current price based off calculation
        S_cur = S_adjust * exp(sqrt(volatility*volatility*time)*gauss_bm);
        //increasing the payoff total by what the max of the price is at that time
        payoff_sum += max(strike_price - S_cur, 0.0);
    }
    //return the average of the payoffs in effect coupling with the risk free return metric
    return (payoff_sum / static_cast<double>(num_sims)) * exp(-risk_free*time);
}



//method used to approximate the cdf for the normal distribution
double normal_cdf(double x) 
{
    //declaring a holder value based off formula
    double k = 1.0 / (1.0 + 0.2316419*x);
    //holder for the usms of the values
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    //checking condition of if it is greater than 0 uses complment if it isn't
    if (x >= 0.0) 
    {
        //returns based off of conditional above using formula for cdf
        return (1.0 - (1.0 / (pow(2 * M_PI, 0.5)))*exp(-0.5*x*x) * k_sum);
    }
    //hits else otherwise
    else 
    {
        //finds the complement then 
        return 1.0 - normal_cdf(-x);
    }
}


//used to find d1 and d2 for the analytical solution of pricing european options
double d_j(int which_one, double underlying_price, double strike_price, double risk_free, double volatility, double time) 
{
    //returns the basic equation that is outlined in the paper for whta these values should be set to based off of which value is passed in for the first parameter
    return (log(underlying_price / strike_price) + (risk_free + (pow(-1, which_one - 1))*0.5*volatility*volatility)*time) / (volatility*(pow(time, 0.5)));
}


//calculates the call price for a european option
double call_price(double underlying_price, double strike_price, double risk_free, double volatility, double time)
{
    //returns the evaluation price for a call based off of formula outlined in paper
    return underlying_price * normal_cdf(d_j(1, underlying_price, strike_price, risk_free, volatility, time)) - strike_price*exp(-risk_free*time) * normal_cdf(d_j(2, underlying_price, strike_price, risk_free, volatility, time));
}


//calculates the put price for a european option
double put_price(double underlying_price, double strike_price, double risk_free, double volatility, double time) 
{
    //returns the evaluation price for a put based off of formula outlined in paper
    return -underlying_price*normal_cdf(-d_j(1, underlying_price, strike_price, risk_free, volatility, time)) + strike_price*exp(-risk_free*time) * normal_cdf(-d_j(2, underlying_price, strike_price, risk_free, volatility, time));
}

//struct that was created in order to store the data surronding a price point in time
struct pricedata 
{
    //holds time
    double time;
    //underlying price
    double underlying_price;
    //strike price
    double strikePrice;
    //risk free rate
    double riskFreeRate;
    //volatility
    double volatility;
    //call price
    double callPrice;
    //put price
    double putPrice;

};

int main(int argc, char* argv[]) 
{
    //these were outstreaks for outputting to a file so it wasn't writting everything to console, worked when testing/debugging
    ofstream out;
    //opening outstream
    out.open("data.txt");

    //MPI::Init(argc, argv);
    // starting the clock to time the program
    double start = MPI::Wtime();
    //setting the number of simulated paths for the monte carlo
    int num_sims = 1000000; 
    //setting underlying price 
    double underlying_price = 100.0;
    //setting the strike price
    double strikePriceStart = 100.0;
    //setting the risk free return rate, 5%
    double r = 0.05; 
    //setting the volatility, this is 20%
    double v = 0.2;    

    //this was implemented but not used really
    //a vector of prices so the simulations could be run at various current asset prices so the trader would
    //have more information as the market shifted
    //would allow for realtime decisions
    //but was moved away from after the serial time was so long
    vector<double> pricingSpan;
    //setting the span for pricing using equal spacing
    pricingSpan = equalSpacing(underlying_price, .25, 1);

    //starting vector to store the pricing data for the monte method
    vector<pricedata> priceVecMonte;
    //starting vectors to store the pricing data for the sto method
    vector<pricedata> priceVecSto;

    //iterate over a span of prices of asset
    for (int j = 0; j < pricingSpan.size(); j++) 
    {
        MPI::Init(argc, argv);
        //this sets the number of times that will be simualted to in this for loop and the simulations
        int timeSpanNum = 1000000;
        //holder for rank
        int rank;
        //holder for num cors
        int numcores;
        //tag didn't do anything
        int tag = 123;
        //a div for knowning how many each core needs to take
        int div;
        //handles all teh remaining datapoints
        int rem;
        //a vector to hold the timespan
        vector <double> timeSpan;
        //populating the timespan using a linspac
        timeSpan = linspace(0.0001, 1, timeSpanNum);

        //declaring a holder for the points that will be sent everywhere
        vector<double>sendpointsMaster;
        //getting each cores rank
        rank = MPI::COMM_WORLD.Get_rank();
        //getting teh size of the comworld
        numcores = MPI::COMM_WORLD.Get_size();

        //here is where the amount that each core gets is assigned
        //based off the number of points over the number of cores
        div = int((timeSpanNum) / numcores);

        //after this the remainder is calculated using mod
        rem = (timeSpanNum) % numcores;

        //here is where we divide up the timespan for each core
        //have master handle this
        if (rank == 0) 
        {
            //iterating over all the cores
            for (int j = 0; j<numcores; j++)
            {
                //creating a vector to store the points that each core must be sent, they all must be the same size when using scatter
                vector<double> sendpoints = vector<double>(div + rem, 0.0001);
                //creating an index for this so the indexing is proper in the timespan and sendpoints
                int index = 0;
                //for loop that assigns values from the timespan
                for (int n = j * div; (n <= (j + 1) * div) && (index < div + rem); n++) 
                {
                    //assignment happening here
                    sendpoints[index++] = timeSpan[n];
                }
                //inserting the generated vector into the end of the master vector that will be scattered
                sendpointsMaster.insert(sendpointsMaster.end(), sendpoints.begin(), sendpoints.end());
            }
        }
        //creating a vector to catch points from our scatter
        vector<double> localPts(div + rem, 0.0);
        //scatter the points to all cores
        MPI::COMM_WORLD.Scatter(&sendpointsMaster[0], div + rem, MPI_DOUBLE, &localPts[0], div + rem, MPI_DOUBLE, 0);
        //this gets rid of any poins that are not actually assigned
        while (!localPts.empty() && localPts.back() <= .00000001) 
        {
            //popping none needed elements off the back
            localPts.pop_back();
        }



        //go over each of the time points
        //each core does this
        for (int i = 0; i < localPts.size(); i++) 
        {
            //creates a holder for monte
            pricedata holderMonte;
            //creates a holder for sto
            pricedata holderSto;
            

            //calls the monte pricing methods here with our parameters
            //first call
            double callMonte = monte_call(num_sims, pricingSpan[j], strikePriceStart, r, v, localPts[i]);
            //then put
            double putMonte = monte_put(num_sims, pricingSpan[j], strikePriceStart, r, v, localPts[i]);
            //all of this is assigning the information into our holder struct
            holderMonte.underlying_price = pricingSpan[j];
            holderMonte.riskFreeRate = r;
            holderMonte.strikePrice = strikePriceStart;
            holderMonte.time = localPts[i];
            holderMonte.volatility = v;
            holderMonte.callPrice = callMonte;
            holderMonte.putPrice = putMonte;

            //after this is assigned then it is pushed back onto the holder that each core has for its pricing data
            priceVecMonte.push_back(holderMonte);

            //calls the sto pricing methods here with our parameters
            //first call
            double callSto = call_price(pricingSpan[j], strikePriceStart, r, v, localPts[i]);
            //then put
            double putSto = put_price(pricingSpan[j], strikePriceStart, r, v, localPts[i]);
            //all of this is assigning the information into our holder struct
            holderSto.underlying_price = pricingSpan[j];
            holderSto.riskFreeRate = r;
            holderSto.strikePrice = strikePriceStart;
            holderSto.time = localPts[i];
            holderSto.volatility = v;
            holderSto.callPrice = callSto;
            holderSto.putPrice = putSto;

            //after this is assigned then it is pushed back onto the holder that each core has for its pricing data
            priceVecSto.push_back(holderSto);


            //this here was outputting evey single pricing point to the console to check and validate that the program was working
            // Finally we output the parameters and prices
            //  cout << "Time:            " << localPts[i] << endl;
            //    cout << "Underlying:      " << pricingSpan[j] << endl;
            //    cout << "Strike:          " << strikePriceStart << endl;
            //    cout << "Risk-Free Rate:  " << risk_free << endl;
            //    cout << "Volatility:      " << volatility << endl;

            //  out << "Monte Carlo Call Price:      " << callMonte << endl;
            //    out << "Monte Carlo Put Price:       " << putMonte << endl;
            //    out << "Sto Call Price:      " << callSto << endl;
            //    out << "Sto Put Price:       " << putSto << endl;
        }

        //declaring a holder for the max sto call price, and assigning it
        double stoCallMax = priceVecSto[0].callPrice;

        //declaring a holder for the min sto call price, and assigning it
        double stoCallMin = priceVecSto[0].callPrice;

        //declaring a holder for the max sto put price, and assigning it
        double stoPutMax = priceVecSto[0].putPrice;

        //declaring a holder for the min sto put price, and assigning it
        double stoPutMin = priceVecSto[0].putPrice;

        //declaring a holder for the max monte call price, and assigning it
        double monteCallMax = priceVecMonte[0].callPrice;
        //declaring a holder for the min monte call price, and assigning it
        double monteCallMin = priceVecMonte[0].callPrice;

        //declaring a holder for the max monte put price, and assigning it
        double montePutMax = priceVecMonte[0].putPrice;
        //declaring a holder for the min monte put price, and assigning it
        double montePutMin = priceVecMonte[0].putPrice;



        //here we find the mins and maxs for each timespan for the sto method
        for (int i = 0; i < priceVecSto.size(); i++) 
        {
            //checking logic for each position in the holder vector and to see if the holders for min/max call/put need to be reassigned
            //checking stoCallMin
            if (priceVecSto[i].callPrice < stoCallMin) 
            {
                //reassignment
                stoCallMin = priceVecSto[i].callPrice;
            }
            //checking stoCallMax
            if (priceVecSto[i].callPrice > stoCallMax) 
            {
                //reassignment
                stoCallMax = priceVecSto[i].callPrice;
            }
            //checking stoPutMin
            if (priceVecSto[i].putPrice < stoPutMin) 
            {
                //reassignment
                stoPutMin = priceVecSto[i].putPrice;
            }
            //checking stoPutMax
            if (priceVecSto[i].putPrice > stoPutMax) 
            {
                //reassignment
                stoPutMax = priceVecSto[i].putPrice;
            }
        }
        //here we find the mins and maxs for each timespan for the monte method
        for (int i = 0; i < priceVecMonte.size(); i++) 
        {
            //checking logic for each position in the holder vector and to see if the holders for min/max call/put need to be reassigned
            //checking monteCallMin
            if (priceVecMonte[i].callPrice < monteCallMin) 
            {
                //reassignment
                monteCallMin = priceVecMonte[i].callPrice;
            }
            //checking monteCallMax
            if (priceVecMonte[i].callPrice > monteCallMax) 
            {
                //reassignment
                monteCallMax = priceVecMonte[i].callPrice;
            }
            //checking montePutMin
            if (priceVecMonte[i].putPrice < montePutMin) 
            {
                //reassignment
                montePutMin = priceVecMonte[i].putPrice;
            }
            //checking montePutMax
            if (priceVecMonte[i].putPrice > montePutMax) 
            {
                //reassignment
                montePutMax = priceVecMonte[i].putPrice;
            }
        }
        //here will gather these from the cores and put them in an ret_array and find the min and max for each method from all of them

        //declaring a vector to send the the min from the sto call
        vector<double> sendMinCallSto;
        //pushing back value
        sendMinCallSto.push_back(stoCallMin);

        //declaring a vector to send the the max from the sto call
        vector<double> sendMaxCallSto;
        //pushing back value
        sendMaxCallSto.push_back(stoCallMax);

        //declaring a vector to send the the min from the sto put
        vector<double> sendMinPutSto;
        //pushing back value
        sendMinPutSto.push_back(stoPutMin);

        //declaring a vector to send the the max from the sto put
        vector<double> sendMaxPutSto;
        //pushing back value
        sendMaxPutSto.push_back(stoPutMax);

        //declaring a vector to send the the min from the monte call
        vector<double> sendMinCallMonte;
        //pushing back value
        sendMinCallMonte.push_back(monteCallMin);

        //declaring a vector to send the the max from the monte call
        vector<double> sendMaxCallMonte;
        //pushing back value
        sendMaxCallMonte.push_back(monteCallMax);

        //declaring a vector to send the the min from the monte put
        vector<double> sendMinPutMonte;
        //pushing back value
        sendMinPutMonte.push_back(montePutMin);

        //declaring a vector to send the the max from the monte put
        vector<double> sendMaxPutMonte;
        //pushing back value
        sendMaxPutMonte.push_back(montePutMax);


        //declaring a global holder
        vector<double>GlobalMinCallSto(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMinCallSto[0], 1, MPI_DOUBLE, &GlobalMinCallSto[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMaxCallSto(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMaxCallSto[0], 1, MPI_DOUBLE, &GlobalMaxCallSto[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMinPutSto(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMinPutSto[0], 1, MPI_DOUBLE, &GlobalMinPutSto[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMaxPutSto(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMaxPutSto[0], 1, MPI_DOUBLE, &GlobalMaxPutSto[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMinCallMonte(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMinCallMonte[0], 1, MPI_DOUBLE, &GlobalMinCallMonte[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMaxCallMonte(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMaxCallMonte[0], 1, MPI_DOUBLE, &GlobalMaxCallMonte[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMinPutMonte(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMinPutMonte[0], 1, MPI_DOUBLE, &GlobalMinPutMonte[0], 1, MPI_DOUBLE, 0);

        //declaring a global holder
        vector<double>GlobalMaxPutMonte(numcores, 0.0);
        //gather from all the cores
        MPI::COMM_WORLD.Gather(&sendMaxPutMonte[0], 1, MPI_DOUBLE, &GlobalMaxPutMonte[0], 1, MPI_DOUBLE, 0);


        //this was outputting to the core in a previous serial version
        /* out << "Sto Call Min: " << stoCallMin<< endl;
        out << "Sto Call Max: " << stoCallMax<< endl;
        out << "Sto Put Min: " << stoPutMin << endl;
        out << "Sto Put Max: " << stoPutMax << endl;

        out << "Monte Call Min: " << monteCallMin << endl;
        out << "Monte Call Max: " << monteCallMax << endl;
        out << "Monte Put Min: " << montePutMin << endl;
        out << "Monte Put Max: " << montePutMax << endl; */


        if (rank == 0) 
        {
            //here is where the master goes thourgh and checks
            //to find the min and max for each option

            //assigning value for a holder
            double globalMinCallSto = GlobalMinCallSto[0];
            //assigning value for a holder
            double globalMaxCallSto = GlobalMaxCallSto[0];
            //assigning value for a holder
            double globalMinPutSto = GlobalMinPutSto[0];
            //assigning value for a holder
            double globalMaxPutSto = GlobalMaxPutSto[0];

            //assigning value for a holder
            double globalMinCallMonte = GlobalMinCallMonte[0];
            //assigning value for a holder
            double globalMaxCallMonte = GlobalMaxCallMonte[0];
            //assigning value for a holder
            double globalMinPutMonte = GlobalMinPutMonte[0];
            //assigning value for a holder
            double globalMaxPutMonte = GlobalMaxPutMonte[0];

            //basically these  iterate over all the items gathered from the cores and they check the values to see if mins and
            //maxs need to be reassigned based off of what was gathered
            for (int i = 0; i < GlobalMinCallSto.size(); i++) 
            {
                //checking condition for min pn sto calls
                if (GlobalMinCallSto[i]<globalMinCallSto)
                {
                    //reassignment
                    globalMinCallSto = GlobalMinCallSto[i];
                }
                //checking condition for max on sto calls
                if (GlobalMaxCallSto[i]>globalMaxCallSto) 
                {
                    //reassignment
                    globalMaxCallSto = GlobalMaxCallSto[i];
                }
                //checking condition for min on sto puts
                if (GlobalMinPutSto[i]<globalMinPutSto) 
                {
                    //reassignment
                    globalMinPutSto = GlobalMinPutSto[i];
                }
                //checking condition for max on sto puts
                if (GlobalMaxPutSto[i]>globalMaxPutSto) 
                {
                    //reassignment
                    globalMaxPutSto = GlobalMaxPutSto[i];
                }
            }
            //repreating the same for monte method
            for (int i = 0; i < GlobalMinCallMonte.size(); i++)
            {
                //checking condition for monte min call
                if (GlobalMinCallMonte[i]<globalMinCallMonte)
                {
                    //reassignment
                    globalMinCallMonte = GlobalMinCallMonte[i];
                }

                //checking condition for monte max call
                if (GlobalMaxCallMonte[i]>globalMaxCallMonte) 
                {
                    //reassignment
                    globalMaxCallMonte = GlobalMaxCallMonte[i];
                }

                //checking condition for monte min put
                if (GlobalMinPutMonte[i]<globalMinPutMonte) 
                {
                    //reassignment
                    globalMinPutMonte = GlobalMinPutMonte[i];
                }

                //checking condition for monte max put
                if (GlobalMaxPutMonte[i]>globalMaxPutMonte) 
                {
                    //reassignment
                    globalMaxPutMonte = GlobalMaxPutMonte[i];
                }
            }

            //after all the values are gone through then it is simply a matter of outputting the information and
            //letting the trader go to work with it and make their specualtive trades
            cout << "The min value of a call using a stochastic method is: " << globalMinCallSto << endl;
            cout << "The max value of a call using a stochastic method is: " << globalMaxCallSto << endl;
            cout << "The min value of a put using a stochastic method is: " << globalMinPutSto << endl;
            cout << "The max value of a put using a stochastic method is: " << globalMaxPutSto << endl;

            cout << "The min value of a call using the monte carlo method is: " << globalMinCallMonte << endl;
            cout << "The max value of a call using the monte carlo method is: " << globalMaxCallMonte << endl;
            cout << "The min value of a put using the monte carlo method is: " << globalMinPutMonte << endl;
            cout << "The max value of a put using the monte carlo method is: " << globalMaxPutMonte << endl;


        }

    }

    //using wtime to end
    double end = MPI::Wtime();
    //calculating the time it took to run
    double elapsed_secs = double(end - start);
    //outputting time
    cout << "Seconds:     " << elapsed_secs << endl;
    //closing output
    out.close();
    //finalize mpi
    MPI::Finalize();
    //done
    return 0;
}
