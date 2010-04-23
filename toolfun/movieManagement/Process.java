// 
// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Process {
	// Constructor
	protected Process(movieData owner, String name, String dateTime)
	{
		owner_ = owner;
		hasChanged_ = false;
	}
	
	// Get the name of the process
	public String getName() {
		return name_;
	}
	
	// Get the date when the process has been finished
	public String getDateTime() {
		return dateTime_;
	}		
	
	// Get the comment the user could have wrote on the process GUI.
	public String getNote() {
		return note_;
	}
	
	// Return whether the process parameters have been changed. The field
	// hasChanged_ needs to be update anytime the process parameters are
	// modified in the process GUI (hasChanged_ = true) or if the process has
	// has just been run with the new parameters (hasChanged_ = false).
	public boolean hasChanged() {
		return hasChanged_;
	}
	
	// abstract methods section
	
	// Get the help text
	public abstract String getHelp();
	
	// make a sanity check of the process
	public abstract void sanityCheck(boolean full);
	
	// protected field section
	protected const movieData owner_;
	protected boolean hasChanged_;
	protected String dateTime_;	
	protected String note_;
	
	// private fields
	private String name_;
}
