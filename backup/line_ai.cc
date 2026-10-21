#include<lv1.hh>
#include<set>
#include<chrono>

// Minimal thread pool implementation
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>

class ThreadPool
{
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
public:
    ThreadPool(size_t n) : stop(false)
    {
        for (size_t i = 0; i < n; ++i)
        {
            workers.emplace_back([this]
            {
                while (true)
                {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                        return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
        }
    }
    template<class F> void enqueue(F&& f)
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.emplace(std::forward<F>(f));
        }
        condition.notify_one();
    }
    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker: workers) worker.join();
    }
};

void lv1::line(string path,uint q0,bool out)
{
    fstream f;
	f.open(path + "/line",fstream::in);
	if(f.is_open())
	{
		f.close();
		return;
	}
	f.close();
	f.open(path + "/line_info",fstream::out | fstream::trunc);
	f<<"t_{f}\nt_{d}\nt_{s}\n#xi_{x}\n#xi_{y}\n";
	f<<"#epsilon_{x}\n#epsilon_{y}\n";
	f<<"#epsilon_{x0}\n#epsilon_{y0}\n";
	f<<"#delta_{0}\nT\n#nu_{s}"<<endl;
	f.close();
	int NL=0, Nc=0, Nr=0, Ncr=0;
	string s_NL = "cat $(find " + path + " -type f -name 'spin.dat' | sed -n '1p') | wc -l";
	string s_Nc = "find " + path + " -type d -name 'c*' | wc -l";
	NL = bash::to_int(s_NL);
	Nc = bash::to_int(s_Nc);
	cout<<"NL = "<<NL<<endl;
	cout<<"Nc = "<<Nc<<endl;

	int nth = 32, q = 0;

    ThreadPool pool(nth);
    std::mutex file_mutex;

    for(int i=0; i<Nc; i++)   
    {
        if(q>=q0)
        {
            pool.enqueue([&, i]()
            {
                array<double,11> ar;
                cout<<"Currently working on c"<<i+1<<endl;
                proton_RF P(path + "/c" + to_string(i+1));
                auto CSx = P.Courant_Snyder(0,out);
                auto CSy = P.Courant_Snyder(1,out);
                auto syn = P.longitudinal(out);
                auto eta = P.phase_slip(out);
                double nus = P.spin_tune(out);

                double kt1 = P.Kt1(), kt2 = P.Kt2(), kt3 = P.Kt3();
                double Cx = P.Cx(), Cy = P.Cy();
                double epx = CSx[0], epy = CSy[0];
                double dta = syn[3], om = syn[4];
                double ps = eta[1];
                ar = {kt1,kt2,kt3,Cx,Cy,epx,epy,dta,om,ps,nus};

                // Write result to file (thread-safe)
                std::lock_guard<std::mutex> lock(file_mutex);
                fstream f2;
                f2.open(path + "/line",fstream::out | fstream::app);
                for(int l=0; l<10; l++) f2<<ar[l]<<"\t";
                f2<<ar[10]<<endl;
                f2.close();
            });
            this_thread::sleep_for(chrono::milliseconds(1)); // Optional: throttle task submission
        }
        else
        {
            cout<<"Thread "<<q%nth<<" not started: skipping line "<<q<<endl;
        }
        q++;
    }
    // ThreadPool destructor waits for all tasks to finish
}