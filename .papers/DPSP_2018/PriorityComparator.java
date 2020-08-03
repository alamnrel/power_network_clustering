import java.util.Comparator;

public class PriorityComparator implements Comparator<double[]>
{    
    @Override
    public int compare(double[] x, double[] y)
    {                
        if (x[2] < y[2])            
            return -1;       
        else if (x[2] > y[2])
            return 1; 
        else 
            return 0;
    }  
}
