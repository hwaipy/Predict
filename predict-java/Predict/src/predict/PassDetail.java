package predict;

import java.util.ArrayList;

/**
 *
 * @author hwaipy
 */
public class PassDetail {
    
    ArrayList<PassPosition> passPositions = new ArrayList<>();

    void append(double daynum, double ele, double azi, double range) {
        passPositions.add(new PassPosition(daynum, ele, azi, range));
    }
}
