#pragma once
/** \brief timer that uses c_clock*/
struct sys_timer{
private:
	std::chrono::high_resolution_clock::time_point t;
	std::chrono::milliseconds ts;
	bool runs;
public:
	sys_timer();	
	void start();	
	void take();
	void stop();
	double get();	
	double stop_and_get();	
	void add_time(std::chrono::milliseconds t);	
	void add_time(sys_timer & timer);
};

sys_timer::sys_timer(){
	ts = std::chrono::milliseconds(0);
	runs = false;
};

void sys_timer::start(){
	t = std::chrono::high_resolution_clock::now();
	runs = true;
};

void sys_timer::take(){
	if(runs){
		std::chrono::high_resolution_clock::time_point  end = std::chrono::high_resolution_clock::now();
		std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - t);
		ts += elapsed;
		t = std::chrono::high_resolution_clock::now();
	}else{
		gedmap_io::print_error("Timer has not been started");
	}
};

void sys_timer::stop(){
	if(runs){
		std::chrono::high_resolution_clock::time_point  end = std::chrono::high_resolution_clock::now();
		std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - t);
		ts += elapsed;
		runs = false;
	}else{
		gedmap_io::print_error("Timer has not been started");
	}
};

double sys_timer::get(){
	uint64_t t =  ts.count();
	return  (double) t / 1000;
};

double sys_timer::stop_and_get(){
	stop();
	return get();
};

void sys_timer::add_time(std::chrono::milliseconds t){
	ts += t;
};

void sys_timer::add_time(sys_timer & timer){
	ts += timer.ts;
};
