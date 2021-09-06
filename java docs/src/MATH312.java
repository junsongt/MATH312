import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Scanner;

public class MATH312 {
    private ArrayList<String> alphabet = new ArrayList<>();

    private String firstFreqLetter = "E";
    private String secondFreqLetter = "T";
    private String thirdFreqLetter = "A";


    public MATH312() {
        alphabet.add("A"); alphabet.add("B"); alphabet.add("C"); alphabet.add("D"); alphabet.add("E");
        alphabet.add("F"); alphabet.add("G"); alphabet.add("H"); alphabet.add("I"); alphabet.add("J");
        alphabet.add("K"); alphabet.add("L"); alphabet.add("M"); alphabet.add("N"); alphabet.add("O");
        alphabet.add("P"); alphabet.add("Q"); alphabet.add("R"); alphabet.add("S"); alphabet.add("T");
        alphabet.add("U"); alphabet.add("V"); alphabet.add("W"); alphabet.add("X"); alphabet.add("Y");
        alphabet.add("Z");
    }

    // global def exponential calculation
    public int power(int a, int d) {
        int product = 1;
        for (int i = 1; i <= d; i++) {
            product = product * a;
        }
        return product;
    }



//===============================================================================================
    //Basics

    // compute GCD
    public int gcd(int a, int b) {
        while (b != 0) {
            int r = a % b;
            a = b;
            b = r;
        }
        return a;
    }

    // compute # of steps taken to get GCD by Euclidean Algorithm
    public int gcdSteps(int a, int b) {
        int step = 0;
        while (b != 0) {
            int r = a % b;
            a = b;
            b = r;
            step = step + 1;
        }
        return step;
    }


    // compute the highest power of base b that n contains
    // n must be positive integer
    public int highestPower(int n, int b) {
        int degree = 0;
        while (n >= power(b, degree)) {
            degree = degree + 1;
        }
        return (degree - 1);
    }


    // compute number of digits of an integer
    public int digitNum(int a) {
//        Integer i = a;
//        String num = i.toString();
//        return num.length();
        return highestPower(a, 10) + 1;
    }


    // list out the binary powers that a given positive integer is composed of;
    // i.e. 100 = 2^6 + 2^5 + 2^2 => [6, 5, 2]
    // i.e. 1109 = 2^10 + 2^6 + 2^4 + 2^2 + 2^0 => [10, 6, 4, 2, 0]
    public ArrayList<Integer> binaryPowerRep(int n) {
        ArrayList<Integer> powerList = new ArrayList<>();
        while (n > 0) {
            int power = highestPower(n, 2);
            powerList.add(power);
            n = n - power(2, power);
        }
        return powerList;
    }


    // compute the modular exponetiation by repeated squares while reducing by modulo
    public int modExp(int a, int e, int n) {
        int result = 1;
        // trivial case: e = 0 => a^0 = 1
        if (e == 0) {
            return (result);
        }

        // when e != 0, then e could be expressed in binary power representation
        // first, list out the binary powers of e
        // then sort the powers from min to max
        ArrayList<Integer> powerList = binaryPowerRep(e);
        powerList.sort(Comparator.naturalOrder());

        int startPower = 0;
        int temp = 1;
        for (int i = 0; i < powerList.size(); i++) {
            int power = powerList.get(i);
            for (int j = startPower; j <= power; j++) {
                if (j == 0) {
                    temp = a;
                } else {
                    temp = power(temp, 2) % n;
                }
            }
            startPower = power + 1;
            result = (result * temp) % n;
        }
        return result;
    }


    // list out solutions 1 <= x <= 50 to the discrete log congruence: b^x = r (mod n)
    public ArrayList<Integer> discreteLog(int b, int r, int n) {
        ArrayList<Integer> logs = new ArrayList<>();
        for (int i = 1; i <= 50; i++) {
            if ((modExp(b, i, n) - r) % n == 0) {
                logs.add(i);
            }
        }
        return logs;
    }


    // check if congruent modulo n
    public boolean cong(int a, int b, int n) {
        return (((a - b) % n) == 0);
    }


    // least non-negative residue
    public int residue(int a, int n) {
        int r = a % n;
        if (r >= 0) {
            return r;
        } else {
            return (n + r);
        }
    }

    // compute Bezout's coefficients, i.e. 6x + 9y = 12 ==> [-4, 4]
    public ArrayList<Integer> bezoutCoeff(int a, int b) {
        int g = gcd(a, b);
        a = a / g;
        b = b / g;
        ArrayList<Integer> v1 = new ArrayList<>();
        ArrayList<Integer> v2 = new ArrayList<>();
        v1.add(1);
        v1.add(0);
        v2.add(0);
        v2.add(1);
        while (b != 1) {
            int r = a % b;
            int q = a / b;
            ArrayList<Integer> v3 = new ArrayList<>();
            v3.add(v1.get(0)- q * v2.get(0));
            v3.add(v1.get(1)- q * v2.get(1));
            v1 = v2;
            v2 = v3;
            a = b;
            b = r;
        }
        return v2;
    }


    // multiplicative inverse
    // compute by extracting the first element of Bezout's coeffs, and make it least non-negative
    public int inverse(int a, int n) {
        if (gcd(a, n) != 1) {
            throw new ArithmeticException();
        }
        ArrayList<Integer> coeff = bezoutCoeff(a, n);
        int inverse = coeff.get(0);
        if (inverse < 0) {
            inverse = n + inverse;
        }
        return inverse;
    }




//================================================================================================
    // ISBN related
    public boolean validISBN(String isbn) {
        int sum = 0;
        for (int i = 0; i < isbn.length(); i++) {
            String s = isbn.substring(i, i+1);
            if (s.equals("X")) {
                s = "10";
            }
            int digit = Integer.parseInt(s);
            sum = sum + digit * (i+1);
        }
        return (sum % 11 == 0);
    }



    public int digitValue(String s) {
        if (s.equals("X")) {
            s = "10";
        }
        int value = Integer.parseInt(s);
        return value;
    }


    // recover corrupted ISBN code with given corrupted position
    public int recoverISBN(String isbn, int pos) {
        int sum = 0;
        for (int i = 0; i < pos -1; i++) {
            String s = isbn.substring(i, i+1);
            int digit = digitValue(s);
            sum = sum + (i + 1) * digit;
        }
        for (int j = pos-1; j < 9; j++) {
            String s = isbn.substring(j, j+1);
            int digit = digitValue(s);
            sum = sum + (j + 2) * digit;
        }
        int a = - (sum * inverse(pos, 11)) % 11;
        int missing = residue(a, 11);

        return missing;
    }


    // recover corrupted ISBN code with corrupted position marked "?"
    public int recoverISBN(String isbn) {
        int sum = 0;
        int length = isbn.length();
        int missingPos = -1;
        for (int i = 0; i < length; i++) {
            String s = isbn.substring(i, i+1);
            if (s.equals("?")) {
                sum = sum + (i + 1) * 0;
                missingPos = i + 1;
            } else {
                int digit = digitValue(s);
                sum = sum + (i + 1) * digit;
            }
        }
        int a = - (sum * inverse(missingPos, 11)) % 11;
        int missing = residue(a, 11);
        return missing;
    }




//================================================================================================
    //Cryptography

    //shift cipher
    // encode with key
    public String shiftEncode(String p, int key) {
        String c = "";
        for (int i = 0; i < p.length(); i++) {
            String sp = p.substring(i, i+1);
            int indexP = alphabet.indexOf(sp);
            int indexC = residue((indexP + key) % 26, 26);
            String sc = alphabet.get(indexC);
            c = c + sc;
        }
        return c;
    }

    // decode with key
    public String shiftDecode(String c, int key) {
        String p = "";
        for (int i = 0; i < c.length(); i++) {
            String sc = c.substring(i, i+1);
            int indexC = alphabet.indexOf(sc);
            int indexP = residue((indexC - key) % 26, 26);
            String sp = alphabet.get(indexP);
            p = p + sp;
        }
        return p;
    }



    // Frequency Analysis
    // give a frequency stats for all letters in a given ciphertext
    public HashMap<String, Integer> letterStats(String s) {
        HashMap<String, Integer> letterMap = new HashMap<String, Integer>();
        for (int i = 0; i < s.length(); i++) {
            String curr = s.substring(i, i+1);
            if (letterMap.containsKey(curr)) {
                int counts = letterMap.get(curr);
                counts = counts + 1;
                letterMap.put(curr, counts);
            } else {
                letterMap.put(curr, 1);
            }
        }
        return letterMap;
    }


    // find the most frequent letter in ciphertext from a letter frequency stats(HashMap)
    public String maxFreqLetter(HashMap<String, Integer> letterMap) {
        String maxletter = "";
        int maxCount = 0;
        for (String key : letterMap.keySet()) {
            int count = letterMap.get(key);
            if (count >= maxCount) {
                maxCount = count;
                maxletter = key;
            }
        }
        return maxletter;
    }

    // list out the most 3 frequent letters in the ciphertext
    public ArrayList<String> mostFrequentLetters(String s) {
        HashMap<String, Integer> letterMap = letterStats(s);
        ArrayList<String> letterRanking = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            String first = maxFreqLetter(letterMap);
            letterRanking.add(first);
            letterMap.remove(first);
        }
        return letterRanking;

//        HashMap<String, Integer> letterMap = new HashMap<String, Integer>();
//        for (int i = 0; i < s.length(); i++) {
//            String curr = s.substring(i, i+1);
//            if (letterMap.containsKey(curr)) {
//                int counts = letterMap.get(curr);
//                counts = counts + 1;
//                letterMap.put(curr, counts);
//            } else {
//                letterMap.put(curr, 1);
//            }
//        }
//        String maxletter = "";
//        int maxCount = 0;
//        for (String key : letterMap.keySet()) {
//            int count = letterMap.get(key);
//            if (count >= maxCount) {
//                maxCount = count;
//                maxletter = key;
//            }
//        }
//        return maxletter;
    }


    // break shift code by frequency analysis to guess key
    public String shiftCrack(String c) {
        String mostFreq = mostFrequentLetters(c).get(0);
        int key = alphabet.indexOf(mostFreq) - alphabet.indexOf(firstFreqLetter);
        return shiftDecode(c, key);
    }



    // affine cipher
    // (a, b) is key
    // encode with enciphering key
    public String affineEncode(String p, int a, int b) {
        String c = "";
        for (int i = 0; i < p.length(); i++) {
            String sp = p.substring(i, i+1);
            int indexP = alphabet.indexOf(sp);
            int indexC = residue((a * indexP + b) % 26, 26);
            String sc = alphabet.get(indexC);
            c = c + sc;
        }
        return c;
    }

    // decode with enciphering key
    public String affineDecode(String c, int a, int b) {
        a = inverse(a, 26);
        b = residue(-b * a, 26);
        String p = "";
        for (int i = 0; i < c.length(); i++) {
            String sc = c.substring(i, i+1);
            int indexC = alphabet.indexOf(sc);
            int indexP = residue((a * indexC + b) % 26, 26);
            String sp = alphabet.get(indexP);
            p = p + sp;
        }
        return p;
    }

    // break affine code by frequency analysis to guess key
    public String affineCrack(String c) {
        ArrayList<String> letterStats = mostFrequentLetters(c);
        String first = letterStats.get(0);
        String second = letterStats.get(1);
        int c1 = alphabet.indexOf(first);
        int c2 = alphabet.indexOf(second);
        int p1 = alphabet.indexOf(firstFreqLetter);
        int p2 = alphabet.indexOf(secondFreqLetter);
        int a = (inverse(residue(p1-p2, 26), 26) * residue(c1-c2, 26)) % 26;
        int b = residue((c1 - a * p1) % 26, 26);
        try {
            return affineDecode(c, a, b);
        } catch(ArithmeticException e) {
            second = letterStats.get(2);
            c2 = alphabet.indexOf(second);
            a = (inverse(residue(p1-p2, 26), 26) * residue(c1-c2, 26)) % 26;
            b = residue((c1 - a * p1) % 26, 26);
            return affineDecode(c, a, b);
        }
    }




    // exponetiation cipher

    // letter grouping by given prime p, producing # of letters per group or block
    // i.e. 2525 < p < 252525, grouping = 2
    public int grouping(int p) {
        int length = digitNum(p);
        if (length % 2 == 1) {
            return (length - 1) / 2;
        } else {
            // creating bound = 252525...25 which has same length as given p
            int bound = 0;
            for (int i = 0; i < length; i = i + 2) {
                bound = bound + 25 * power(10, i);
            }
            if (p >= bound) {
                return (length / 2);
            } else {
                return (length / 2 - 1);
            }
        }
    }



    // e is encryptrion key
    public String expEncode(String p, int e, int n) {
        String cipher = "";
        int length = p.length();
        int unitLength = grouping(n);
        int r = length % unitLength;
        // check if the string needs adding additional letters at the end to make group
        for (int k = 0; k < r; k++) {
            p = p + "X";
        }

        for (int i = 0; i < p.length(); i = i + unitLength) {
            //extract the block
            String block = p.substring(i, i+unitLength);
            // temp is each block's numeric value, i.e. "EH" = 0407 = 407
            int temp = 0;
            for (int j = 0; j < unitLength; j++) {
                // extract single letter in each block
                String letter = block.substring(j, j+1);
                int index = alphabet.indexOf(letter);
                // power is the exponent of bases(10^2k) forming numeric value
                // i.e. "EH" = 0407 = 407 = 04 * 10^2 + 07 * 10^0
                int power = (unitLength - 1 - j) * 2;
                // accumulate block numeric value by each letter value
                temp = temp + index * power(10, power);
            }

            Integer result = modExp(temp, e, n);
            String sc = result.toString();
            cipher = cipher + " " + sc;
        }
        return cipher;
    }


//    public String expEncode(String p, int e, int n) {
//        String c = "";
//        for (int i = 0; i < p.length(); i++) {
//            String sp = p.substring(i, i+1);
//            int indexP = alphabet.indexOf(sp);
//            int indexC = modExp(indexP, e, n);
//            String sc = alphabet.get(indexC);
//            c = c + sc;
//        }
//        return c;
//    }


    // e is encryptrion key, decryption key d = e^-1 (mod n-1) by Fermat's little theorem
    public String expDecode(String c, int e, int n) {
        String plaintext = "";
        // find decryption key d
        int d = inverse(e, n-1);
        int unitLength = grouping(n);

        // ciphertext given may contain blank spaces to indicate blocks, i.e. "1203 1175 2307 0103"
        // extract blocks from the ciphertext by cursor stopping at a blank space while deciphering
        int i = 0;
        int end = c.length();
        while (i < end) {
            int cursor = i;
            boolean pause = false;
            while(!pause & cursor < end) {
                String curr = c.substring(cursor, cursor + 1);
                if (curr.equals(" ")) {
                    pause = true;
                } else {
                    cursor = cursor + 1;
                }
            }
            String block = c.substring(i, cursor);
            int temp = Integer.parseInt(block);

            // decipher each block from ciphertext
            // result = direct numeric value of block after applying decryption key d
            int result = modExp(temp, d, n);
            // translate block numeric value to corresponding letters
            for (int u = unitLength; u >= 0; u = u - 2) {
                // power is the exponent of bases(10^2k) forming result
                // 191301 = 19*10^4 + 13*10^2 + 01*10^0
                // 20511 = 020511 = 02*10^4 + 05*10^2 + 11*10^0
                int power = (u - 1) * 2;
                // extract 2 digits from result at a time to be index
                int index = result / power(10, power);
                String sp = alphabet.get(index);
                plaintext = plaintext + sp;
                // prepare for the next 2 digits extraction
                result = result - index * power(10, u);
            }
            // start extracting the next block
            i = cursor + 1;
        }
        return plaintext;

    }

//    public String expDecode(String c, int e, int n) {
//        String p = "";
//        int d = inverse(e, n-1);
//        for (int i = 0; i < p.length(); i++) {
//            String sc = c.substring(i, i+1);
//            int indexC = alphabet.indexOf(sc);
//            int indexP = modExp(indexC, d, n);
//            String sp = alphabet.get(indexP);
//            p = p + sp;
//        }
//        return p;
//    }







//=========================================================================================
    //primitive roots

    // coprime congruence classes
    // reduced residue system
    public ArrayList<Integer> coprimeClass(int n) {
        ArrayList<Integer> classes = new ArrayList<>();
        for (int i = 0; i <= n-1; i++) {
            if (gcd(i, n) == 1) {
                classes.add(i);
            }
        }
        return classes;
    }

    // Euler phi-function
    public int eulerPhi(int n) {
        return coprimeClass(n).size();
    }


    // compute order of integer a (mod n)
    // ord = 0 ====> ord = phi(n)
    public int order(int a, int n) {
        int i = 1;
        while (modExp(a, i, n) != 1) {
            i = i + 1;
        }
        return i;
    }


    // primitive roots (mod n)
    public ArrayList<Integer> primitiveRoot(int n) {
        ArrayList<Integer> classes = coprimeClass(n);
        ArrayList<Integer> roots = new ArrayList<>();
        int phi = eulerPhi(n);
        for (int a : classes) {
            int d = order(a, n);
            if (d == phi) {
                roots.add(a);
            }
        }
        return roots;
    }


    // compute the index i such that root^i = a (mod n) where a is from integers with gcd(a, n) = 1
    // |a| = phi(n)
    public int ind(int root, int a, int n) {
        int index = 0;
        while (modExp(root, index, n) != a) {
            index = index + 1;
        }
        return index;
    }

    // display the index table for every a in {a | gcd(a, n) = 1}
    public HashMap<Integer, Integer> indexTable(int root, int n) {
        ArrayList<Integer> reduceResidues = coprimeClass(n);
        HashMap<Integer, Integer> table = new HashMap<>();
        for (int a : reduceResidues) {
            int index = ind(root, a, n);
            table.put(a, index);
        }
        return table;
    }






//==========================================================================================
    public static void main(String[] args) {
        MATH312 math312 = new MATH312();

        Scanner input = new Scanner(System.in);
        long t1 = System.currentTimeMillis();

//        System.out.println(math312.gcdSteps(23, 17));
//        System.out.println(math312.gcdSteps(15673, 2532));
//        System.out.println(math312.cong(5, 17, 13));
//        System.out.println(math312.cong(5, 18, 13));
//        System.out.println(math312.inverse(7, 6));


//        System.out.println(math312.validISBN("0201065614"));
//        System.out.println(math312.validISBN("354019102X"));
//        System.out.println(math312.validISBN("0198538049"));
//        System.out.println(math312.validISBN("0198738049"));

//        System.out.println(math312.digitNum(1000));
//        System.out.println(math312.digitNum(1352));
//        System.out.println(math312.digitNum(9999));
//        System.out.println(math312.digitNum(10000));
//        System.out.println(math312.digitNum(1));
//        System.out.println(math312.digitNum(9));
//        System.out.println(math312.digitNum(0));

//        System.out.println(math312.recoverISBN("019838049", 5));
//        System.out.println(math312.recoverISBN("915542126", 9));
//        System.out.println(math312.recoverISBN("26105073X", 1));
//        System.out.println(math312.recoverISBN("0198?38049"));
//        System.out.println(math312.recoverISBN("91554212?6"));
//        System.out.println(math312.recoverISBN("?26105073X"));



//        System.out.println(math312.shiftEncode("HELLO", 3));
//        System.out.println(math312.shiftDecode("PHHWDWIRXU", 3));
//        System.out.println(math312.shiftDecode("LFDPHLVDZLFRQTXHUHG", 3));
//        System.out.println(math312.mostFrequentLetters("WSOHQSFSOZSJSFMROM"));
//        System.out.println(math312.mostFrequentLetters("KYVMRCLVFWKYVBVPZJJVMVEKVVE"));
//        System.out.println(math312.shiftCrack("WSOHQSFSOZSJSFMROM"));
//        System.out.println(math312.shiftCrack("KYVMRCLVFWKYVBVPZJJVMVEKVVE"));

//        System.out.println(math312.affineDecode("RTOLKTOIK", 3, 24));
//        System.out.println(math312.affineEncode("THERIGHTCHOICE", 15, 14));



//        System.out.println(math312.binaryPowerRep(1109));

//        System.out.println(math312.modExp(504, 1109, 2881));
//        System.out.println(math312.modExp(1874, 1109, 2881));
//        System.out.println(math312.modExp(347, 1109, 2881));
//        System.out.println(math312.modExp(515, 1109, 2881));
//        System.out.println(math312.modExp(2088, 1109, 2881));
//        System.out.println(math312.modExp(2356, 1109, 2881));
//        System.out.println(math312.modExp(736, 1109, 2881));
//        System.out.println(math312.modExp(468, 1109, 2881));
//        System.out.println(math312.discreteLog(20, 24, 29));

//        System.out.println(math312.affineCrack("MJMZKCXUNMGWIRYVCPUWMPRRWGMIOPMSNYSRYRAZPXMCDWPRYEYXD"));
//        System.out.println(math312.affineCrack("PJXFJSWJNXJMRTJFVSUJOOJWFOVAJRWHEOFJRWJODJFFZBJF"));
        // next one is very tricky!!!
//        System.out.println(math312.affineCrack("WEZBFTBBNJTHNBTADZQETGTYRBZAJNANOOZATWGNABOVGFNWZVA");


//        System.out.println(math312.primitiveRoot(15));
//        System.out.println(math312.primitiveRoot(16));
//        System.out.println(math312.primitiveRoot(19));
//        System.out.println(math312.primitiveRoot(27));
//        System.out.println(math312.primitiveRoot(9));

        // 6 and 8128 are perfect numbers
//        System.out.println(math312.primitiveRoot(6));
//        System.out.println(math312.primitiveRoot(8128));

//        System.out.println(math312.letterStats("PJXFJSWJNXJMRTJFVSUJOOJWFOVAJRWHEOFJRWJODJFFZBJF"));
//        System.out.println(math312.indexTable(5, 23));
//        System.out.println(math312.indexTable(7, 22));


//        System.out.println(math312.grouping(1300));
//        System.out.println(math312.grouping(2524));
//
//        System.out.println(math312.grouping(2526));
//        System.out.println(math312.grouping(25251));
//        System.out.println(math312.grouping(252524));
//
//        System.out.println(math312.grouping(252526));
//        System.out.println(math312.grouping(300000));


//        System.out.println(math312.expEncode("DEEPYOGURT", 11, 2621));
//        System.out.println(math312.expDecode("65 415 1323 1567 150", 11, 2621));







        long t2 = System.currentTimeMillis();
        System.out.println(t2-t1 + "ms");
    }

}
