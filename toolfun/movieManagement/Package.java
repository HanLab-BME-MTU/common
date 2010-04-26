// The first thing to do when a package control panel is opened is to check
// whether an object of the same class exists in the list of packages.
//
// * If the package already exists:
// - a pop up windows should notify the user that the previous package information
//   wil be used.
//
// * otherwise:
// - Call the constructor of the pacakge. Example:
//
//   newPackage = BiosensorPackage(movieData);
//
// In both case, call the sanity check. Any exception thrown during that check
// set the flag icons (not run, problem, dependency problem) next to each process.
//
//   try {
//     newPackage.sanityCheck(full);
//   } catch(ExceptionProcesses) {
//     // display warning icons according to exceptions thrown
//   } catch(otherTypeOfExceptions) {
//     // ...
//   } finally {
//     // No exception caught: every process has to have a "green check" icon
//   }
//	
// - If no critical error occured and the package was new, add the package to the
// movieData package list:
//
//   movieData.addPackage(newPackage);
//
// This class defines the abstract class Package from which every user-defined
// package will inherit.

public abstract class Package {
	// Constructor
	protected Package(movieData owner, String name)
	{
		owner_ = owner;
		name_ = name;
		
		// initialize the process dependency matrix.
		depMatrix_ = getProcessDependenyMatrix_();
		
		// initialize localProcesses_ with the list of required process. Each
		// process are built with default parameters.
		localProcesses_ = getDefaultProcessList_();
		
		// the number of processes required by the package
		int nLocalProcesses = length(localProcesses_);
		
		// the number of total processes stored in the movieData
		int nGlobalProcesses = length(owner_.processes_);
		
		// Check whether a required process has not been already computed by
		// another Package. If so, a pop window should ask the user whether he
		// wants to use the already computed process or he prefers to create a
		// new one. This pop up window should have 2 buttons 'Choose...'
		// and 'Create'. 'Choose' will open another window with the list of
		// existing related processes. 'Create' will create a new process with
		// default parameter.
		
		boolean[] isSame = new int[nProcesses];
		
		// Iterate over every required process and check whether processes of
		// the same class exists
		for (int i = 0; i < nLocalProcesses; ++i) {
			
			String className = className(localProcesses_[i]);
			
			for (int j = 0; j < nGlobalProcesses; ++j) {
				isSame = classname(owner_.processes_[j]) == className;
			}
			
			if (any(isSame)) {
				// At least one process is of the same class as the required ith
				// process. This is when we need to show the popup window with
				// "choose..." and "create" button (see header at the top of the
				// file). Every localProcesses_[j] where isSame[j] == 1 needs to
				// be displayed as a potential choice. 
				
				// TODO: ask user
				
				// test the user answer
				if (response == "choose..." && chosenProcess ~= end) {
					// Replace the default process by the one chosen by user
					localProcesses_[i] = owner_.processes_[chosenProcess];
				} else {
					// Add the new process into owner's process list. 
					owner_.processes_.add(localProcesses_[i]);
				}
			}
		}		
	}
	
	// Get the name of the package
	public String getName() {
		return name_;
	}
	
	// make a sanity check of the package
	public abstract void sanityCheck(boolean full);
	
	// get the list of processes required by the package. These processes contains
	// the default parameters to be used in the algorithm.
	protected abstract Process[] getDefaultProcessList_();
	
	// This method return the matrix of process dependency. It represents the
	// directed graph of processes. Only the lower triangular part of the matrix
	// is considered. For every i > j, mat(i, j) == 1 means the ith process
	// depends on the jth process. If there is no dependency between processes,
	// return a zero matrix, where the size is equal to the number of processes.
	// Any process should have selft dependency, i.e. there should not be any
	// mat(i,i) == 1. The type of the matrix should be logical, i.e.
	// mat = false(nProcess); mat(i,j) = true, etc.
	protected abstract Matrix getProcessDependenyMatrix_();
	
	// This method traverses the list of local processes and call for each process
	// its sanityCheck method. The list of processes is traversed using a depth
	// first search (DFS) algorithm to check dependency between processes.
	protected void checkProcesses_(boolean full) {
		
		// First, call sanityCheck method on each process and store any exception
		// in an array of exception.
		
		int nProcesses = length(localProcesses_); 
		
		// For each process, one needs to store every possible reasons a process
		// is wrong and display them onto the icon (ok, warning triangle, wrong)
		// next to the process.
		//
		// A process status can be alter for 3 reasons:
		// 1. the process itself has a problem
		// 2. the parameters in the process setting panel have changed
		// 3. the process depends on processes which also have a problem.
		//
		// We define an array of vector (Java). In Matlab it will be a cell array,
		// i.e.
		// processExceptions = cell(nProcesses, 1);
		//
		// If a process doesn't have any problem, the cell will be empty.
		
		Vector[] processExceptions = new Vector[nProcesses];
		
		// 1. Check whether the process itself has a problem
		
		for (int i = 0; i < nProcesses; ++i) {
			try {
				localProcesses_[i].sanityCheck(full);
			}
			catch (Exception) {
				processExceptions[i].add(Exception);
			}
		}
		
		// 2. Check whether the process parameters have changed
		
		for (int i = 0; i < nProcesses; ++i) {
			if (localProcesses_[i].hasChanged()) {
				newException = new Exception("Some of the process parameters
											 have changed since previous run.");
				
				processExceptions[i].add(newException);
			}
		}
		
		// 3. Check if a process depends on processes which also have a problem.
		
		boolea[] processVisited = new false[nProcesses];
		
		for (int i = 0; i < nProcesses; ++i) {
			if (!processVisited[i])
				dfs_(i, processExceptions, processVisited);
		}
		
		// Finally, if at least one of the process has a problem, we throw the
		// processExeptions cell array (the way it is done in Matlab is specific
		// so I just write something simple here: 
		if any(processExceptions) {
			throws processExceptions;
		}
	}
	
	private void dfs_(int iProcess, Vector[] processExceptions, processVisited) {
		
		processVisited[iProcess] = true;
		
		for (int j = 0; j < iProcess; ++j) {
			// does process iProcess depend on j ?
			if (depMatrix_[iProcess, j] == 1) {
				// ok, process iProcess depends on process j. Do we need to explore process j?
				if (!processVisited[j]) {
					// yes
					dfs_(j, processExceptions, processVisited);
				}

				// if jth process has a problem, we add an exception to the
				// ith process list of exception since ith process depends
				// on the jth process.
				if (processExceptions[j].length != 0) {
					newException = new Exception("Process " + iProcess +
												 " depends on process " +
												 j + "which has problem.");
					
					processExceptions[iProcess].add(newException);
				}
			}
		}
	}
	
	// private fields
	private String name_;
	
	// protected fields section
	protected movieData owner_;
	// dependency matrix
	protected Matrix depMatrix_;
	// list of required processes
	protected Process[] localProcesses_;
}
