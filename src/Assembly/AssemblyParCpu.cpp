#include "AssemblyParCpu.h"

AssemblyParCpu::AssemblyParCpu(Sparse& s) 
: AssemblyParCpu(s, std::thread::hardware_concurrency()-1)
{}

AssemblyParCpu::AssemblyParCpu(Sparse& s, int numberOfCores) 
: AssemblySingleCpu(s), concurentThreadsSupported(numberOfCores)
{
    Log::Logger().Info("Assembly par CPU created with " + std::to_string(concurentThreadsSupported) + " threads");
    simulation_per_thread();
}

// SubVecSize: number of subvector required
int SubVecSize(int n)
{
    if (n > 1)
        return (std::ceil(static_cast<float>(n/2.0)) + SubVecSize(n/2));
    else
        return 0;
}
void AssemblyParCpu::simulation_per_thread()
{
    int simulationSize = v_dofi.size();
    size_t subVectorSize = concurentThreadsSupported + SubVecSize(concurentThreadsSupported);
    subVector.resize(subVectorSize);
    int loadPerThread =  static_cast<int>(simulationSize/concurentThreadsSupported);
    int loadLastThread = static_cast<int>(simulationSize/concurentThreadsSupported + simulationSize%concurentThreadsSupported);
    for (int i = 0; i < concurentThreadsSupported - 1; i++)
    {
        subVector[i] = Pair(i*loadPerThread,(i+1)*loadPerThread);
    }
    subVector[concurentThreadsSupported-1] = Pair((concurentThreadsSupported-1)*loadPerThread,(concurentThreadsSupported-1)*loadPerThread + loadLastThread);
    
    //
    size_t counter = concurentThreadsSupported;
    size_t max = concurentThreadsSupported;
    size_t min = 1;
    for (size_t i = min; i < max; i = i+2)
    {
        subVector[counter++] = Pair(std::min(subVector[i-1].x, subVector[i].x),
                                    std::max(subVector[i-1].y, subVector[i].y));
    }
    if (max%2) 
        std::rotate(subVector.begin() + max - 1 , subVector.begin() + max  , subVector.begin() + counter);


    for (auto& ii : subVector)
    {
        std::cout << ii << "\n";
    }
    // building subvector for merging 
/*     size_t counter = concurentThreadsSupported;
    size_t max = concurentThreadsSupported;
    size_t min = 1;
    size_t col = 1;
    while (max >= 2) {
        for (size_t i = min; i < max; i = i+2) {
            subVector[counter++] = Pair(subVector[i-1].x,subVector[i].y);
        }
        max = counter - max;
        min = counter;
    }
    if (max%2) {
        Log::Logger().Info(counter);
        subVector[counter] = Pair(subVector[counter-1].x,subVector[counter-2].y);
    }
      */
}

AssemblyParCpu::~AssemblyParCpu() 
{
    Log::Logger().Info("Assembly par CPU destroyed");
}

void AssemblyParCpu::eachSort(int ii)
{
    std::sort(h_indices.begin() + subVector[ii].x, 
              h_indices.begin() + subVector[ii].y, 
              [this](size_t i, size_t j) {
                if (v_dofi[i] == 0 )
                    return false;
                if (v_dofi[j] == 0)
                    return true;
                if (v_dofi[i] == v_dofi[j])
                    return v_dofj[i] < v_dofj[j];
                return v_dofi[i] < v_dofi[j];
                });
}

// fix the merge to run in parallel
void AssemblyParCpu::merge(int ii)
{
    const auto comparison = [this](size_t i, size_t j) {
    if (v_dofi[i] == 0 )
        return false;
    if (v_dofi[j] == 0)
        return true;
    if (v_dofi[i] == v_dofi[j])
        return v_dofj[i] < v_dofj[j];
    return v_dofi[i] < v_dofi[j];
    };
    
    std::merge(h_indices.begin() + subVector[ii].x, h_indices.begin() + subVector[ii].y,
               h_indices.begin() + subVector[ii+1].x, h_indices.begin() + subVector[ii+1].y, 
               h_indices_new.begin() + std::min(subVector[ii].x, subVector[ii+1].x), 
               comparison);
}

void AssemblyParCpu::sort()
{
    Timer timer("Time spend in ParSort -> ");
    std::thread t[concurentThreadsSupported]; 
    for (int i = 0; i < concurentThreadsSupported; i++)
    {
        t[i] = std::thread(&AssemblyParCpu::eachSort,this,i);
    }
    for (int i = 0; i < concurentThreadsSupported; i++)
    {
        t[i].join();
    }
    h_indices_new = h_indices;
    merge(0);
    //subVector[3] = Pair(subVector[0].x,subVector[0+1].y);
    h_indices = h_indices_new;
    merge(2);
    //subVector[4] = Pair(subVector[3].x,subVector[2].y);
    h_indices = h_indices_new;

}
void AssemblyParCpu::calculateAssembly()
{
    Timer *t = new Timer("Time in assembly par cpu -> ");
    //sort();
    //nnzFinder();
    //addDuplicates();
    //eraseDuplicands();
    apply_permutation_in_place(v_value, h_indices);
    apply_permutation_in_place(v_dofi, h_indices);
    apply_permutation_in_place(v_dofj, h_indices);
    delete t;
    //Log::Logger().Info(h_rowPtr.back());
    //stiffMat.valueSize = h_rowPtr.back();
    for (size_t i = 0; i < stiffMat.valueSize; i++) {
        stiffMat.value[i] = v_value[i];
        stiffMat.i[i] = v_dofi[i];
        stiffMat.j[i] = v_dofj[i];
    }
}

// pair class
Pair::Pair(int a, int b)
: x(a), y(b)
{}

Pair::Pair()
{
    x = 0;
    y = 0;
}

// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Pair const& a) 
{
    return out << "<" << a.x << ", " << a.y << ">";
}