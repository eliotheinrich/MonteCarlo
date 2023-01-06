#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <string>
#include <map>

class Sample {
    private:
        double mean;
        double std;
        double num_samples;

    public:
        Sample(double mean);
		Sample(double mean, double std, int num_samples);

        double get_mean();
		void set_mean(double mean);

        double get_std();
		void set_std(double std);

        double get_num_samples();
		void set_num_samples(int num_samples);

        Sample combine(Sample* other);
};

class DataSlide {
	private:
		std::map<std::string, int> data_int;
		std::map<std::string, double> data_double;
		std::map<std::string, std::vector<Sample>> data;

	public:
		DataSlide();
		void add_int(std::string s, int i);
		void add_double(std::string s, double f);
		void add_data(std::string s);
		void push_data(std::string s, Sample sample);
};

class DataFrame {
	private:
		std::vector<DataSlide> slides;


	public:
		DataFrame();
		DataFrame(std::vector<DataSlide> slides);
		void add_slide(DataSlide ds);
		void save(std::string filename);
};


#endif