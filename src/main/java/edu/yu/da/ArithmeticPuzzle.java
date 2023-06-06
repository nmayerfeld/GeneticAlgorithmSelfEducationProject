package edu.yu.da;


import java.util.*;
public class ArithmeticPuzzle extends ArithmeticPuzzleBase{

    public class Solution implements SolutionI{
        private String aug;
        private String add;
        private String sum;
        private List<Character> sol;
        private int generation=0;
        public Solution(String aug, String add, String sum){
            this.aug=aug;
            this.add=add;
            this.sum=sum;
            this.sol=new ArrayList<>();
        }
        public void setSolution(char[] c){
            for(int i=0;i<c.length;i++){
                if(c[i]!='\u0000'){
                    sol.add(c[i]);
                }
                else{
                    sol.add(' ');
                }
            }
        }
        public void nextGeneration(){
            this.generation++;
        }
        @Override
        public List<Character> solution() {
            return  this.sol;
        }
        @Override
        public String getAugend() {
            return this.aug;
        }
        @Override
        public String getAddend() {
            return this.add;
        }
        @Override
        public String getSum() {
            return this.sum;
        }
        @Override
        public int nGenerations() {
            return this.generation;
        }
    }
    public class IndexMaxPQ<Key extends Comparable<Key>> implements Iterable<Integer> { //credit sedgewick and wayne- https://algs4.cs.princeton.edu/24pq/IndexMaxPQ.java
        private int maxN;        // maximum number of elements on PQ
        private int n;           // number of elements on PQ
        private int[] pq;        // binary heap using 1-based indexing
        private int[] qp;        // inverse of pq - qp[pq[i]] = pq[qp[i]] = i
        private Key[] keys;      // keys[i] = priority of i
        public IndexMaxPQ(int maxN) {
            if (maxN < 0) throw new IllegalArgumentException();
            this.maxN = maxN;
            n = 0;
            keys = (Key[]) new Comparable[maxN + 1];    // make this of length maxN??
            pq = new int[maxN + 1];
            qp = new int[maxN + 1];                   // make this of length maxN??
            for (int i = 0; i <= maxN; i++)
                qp[i] = -1;
        }
        public int maxIndex() {
            if (n == 0) throw new NoSuchElementException("Priority queue underflow");
            return pq[1];
        }
        public boolean isEmpty() {
            return n == 0;
        }
        public boolean contains(int i) {
            validateIndex(i);
            return qp[i] != -1;
        }
        public int size() {
            return n;
        }
        public void insert(int i, Key key) {
            validateIndex(i);
            if (contains(i)) throw new IllegalArgumentException("index is already in the priority queue");
            n++;
            qp[i] = n;
            pq[n] = i;
            keys[i] = key;
            swim(n);
        }
        public int delMax() {
            if (n == 0) throw new NoSuchElementException("Priority queue underflow");
            int max = pq[1];
            exch(1, n--);
            sink(1);

            assert pq[n + 1] == max;
            qp[max] = -1;        // delete
            keys[max] = null;    // to help with garbage collection
            pq[n + 1] = -1;        // not needed
            return max;
        }
        public void changeKey(int i, Key key) {
            validateIndex(i);
            if (!contains(i)) throw new NoSuchElementException("index is not in the priority queue");
            keys[i] = key;
            swim(qp[i]);
            sink(qp[i]);
        }

        // throw an IllegalArgumentException if i is an invalid index
        private void validateIndex(int i) {
            if (i < 0) throw new IllegalArgumentException("index is negative: " + i);
            if (i >= maxN) throw new IllegalArgumentException("index >= capacity: " + i);
        }

        private boolean less(int i, int j) {
            return keys[pq[i]].compareTo(keys[pq[j]]) < 0;
        }

        private void exch(int i, int j) {
            int swap = pq[i];
            pq[i] = pq[j];
            pq[j] = swap;
            qp[pq[i]] = i;
            qp[pq[j]] = j;
        }

        private void swim(int k) {
            while (k > 1 && less(k / 2, k)) {
                exch(k, k / 2);
                k = k / 2;
            }
        }

        private void sink(int k) {
            while (2 * k <= n) {
                int j = 2 * k;
                if (j < n && less(j, j + 1)) j++;
                if (!less(k, j)) break;
                exch(k, j);
                k = j;
            }
        }

        public Iterator<Integer> iterator() {
            return new HeapIterator();
        }

        private class HeapIterator implements Iterator<Integer> {
            // create a new pq
            private IndexMaxPQ<Key> copy;

            // add all elements to copy of heap
            // takes linear time since already in heap order so no keys move
            public HeapIterator() {
                copy = new IndexMaxPQ<Key>(pq.length - 1);
                for (int i = 1; i <= n; i++)
                    copy.insert(pq[i], keys[pq[i]]);
            }

            public boolean hasNext() {
                return !copy.isEmpty();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }

            public Integer next() {
                if (!hasNext()) throw new NoSuchElementException();
                return copy.delMax();
            }
        }
    }
    private IndexMaxPQ<Integer> maxFit;
    private char[][] individuals;
    private int[] fitness;
    private String augend;
    private String addend;
    private String sum;
    private int initialPopulationSize;
    private int maxGenerations;
    private double mutationProbability;
    private double crossoverProbability;
    private HashSet<Character> uniqueChars;
    private int solutionIndex;
    private boolean foundSolution=false;
    private int totalFitness=0;

    /**
     * Constructor.  Specifies the arithmetic puzzle to be solved in terms of an
     * augend, addend, and sum.
     * <p>
     * Representation: all characters are in the range A-Z, with each letter
     * represents a unique digit.  The puzzle to be solved specifies that the
     * augend and addend (each representing a number in base 10) sum to the
     * specified sum (also a number in base 10).  Each of these numbers is
     * represented with the most significant letter (digit) in position 0, next
     * most significant letter (digit) in position 1, and so on.  The numbers
     * need not be the same length: an "empty" digit is represented by the
     * "space" character.
     * <p>
     * Addition: Augend + Addend = Sum
     *
     * @param augend
     * @param addend
     * @param sum
     */
    public ArithmeticPuzzle(String augend, String addend, String sum) {
        super(augend, addend, sum);
        this.augend=augend;
        this.addend=addend;
        this.sum=sum;
        this.uniqueChars=new HashSet<>();

        String total=augend+addend+sum;
        for(int i=0;i<total.length();i++){
            char c=total.charAt(i);
            if(!this.uniqueChars.contains(c)){
                this.uniqueChars.add(c);
            }
        }
    }
    @Override
    public SolutionI solveIt(GeneticAlgorithmConfig gac) {
       this.initializeVars(gac);
       this.createRandomPopulation();
       //System.out.println("line 329");
       this.initializeFitnesses();
       Solution s=new Solution(this.augend,this.addend,this.sum);
       if(foundSolution){
           s.setSolution(this.individuals[this.solutionIndex]);
       }
       while(s.nGenerations()<=maxGenerations&&!this.foundSolution){
           //System.out.println("entering loop");
           s.nextGeneration();
           this.processGeneration(gac.getSelectionType(),gac.getMutationProbability(),gac.getCrossoverProbability());
           if(this.foundSolution){
               s.setSolution(this.individuals[this.solutionIndex]);
           }
       }
       return s;
    }
    private void processGeneration(GeneticAlgorithmConfig.SelectionType st, double mutProb,double crossProb){
        Map<Integer,Integer>indexToIndexForCrossover=new HashMap<>();
        if(st.equals("ROULETTE")) indexToIndexForCrossover=this.rouletteSelectionForCrossover();
        else indexToIndexForCrossover=this.tournamentSelectionForCrossover();
        //System.out.println("hit line 348");
        this.performCrossover(indexToIndexForCrossover);
        //System.out.println("hit line 350");
        this.performMutation();
        //this.print();
        //System.out.println("processed generation");
    }
    private void mutate(int index){
        int oldFitness=this.fitness[index];
        char[] c=this.individuals[index];
        /*System.out.println("\n mutating the following array");
        for(char ch:c){
            System.out.print(ch+" ");
        }
        System.out.println();*/
        Set<Integer> indicesSelected=new HashSet<>();
        boolean notFound=true;
        int index1 = 0,index2 =0,index3=0,index4 =0;
        while(notFound){
            index1= (int) (Math.random()*(10));
            if(!indicesSelected.contains(index1)&&c[index1]!='\u0000'){
                notFound=false;
                indicesSelected.add(index1);
            }
        }
        notFound=true;
        while(notFound){
            index2= (int) (Math.random()*(10));
            if(!indicesSelected.contains(index2)){
                notFound=false;
                indicesSelected.add(index2);
            }
        }
        notFound=true;
        while(notFound){
            index3= (int) (Math.random()*(10));
            if(!indicesSelected.contains(index3)&&c[index1]!='\u0000'){
                notFound=false;
                indicesSelected.add(index3);
            }
        }
        notFound=true;
        while(notFound){
            index4= (int) (Math.random()*(10));
            if(!indicesSelected.contains(index4)){
                notFound=false;
                indicesSelected.add(index4);
            }
        }
        char temp=c[index1];
        c[index1]=c[index2];
        c[index2]=temp;
        temp=c[index3];
        c[index3]=c[index4];
        c[index4]=temp;
        this.individuals[index]=c;
        /*System.out.println("printing array post mutation");
        for(char cha: this.individuals[index]){
            System.out.print(cha+" ");
        }
        System.out.println();*/
        int newFitness=this.getFitness(index);
        this.fitness[index]=newFitness;
        this.totalFitness=this.totalFitness-oldFitness+newFitness;
        this.maxFit.changeKey(index,newFitness);
        if(newFitness==0){
            this.foundSolution=true;
            this.solutionIndex=index;
        }
    }
    private void performMutation(){
        List<Integer> selectedForMutation=this.selectForMutation();
        for(Integer i:selectedForMutation){
            double d=Math.random();
            if(d>this.mutationProbability){
                continue;
            }
            //System.out.println("\n mutating at index: "+i);
            this.mutate(i);
        }
    }
    private List<Integer> selectForMutation(){
        List<Integer> list=new ArrayList<>();
        Set<Integer> indicesSelected=new HashSet<>();
        int numForCrossover=this.initialPopulationSize/20;
        if(numForCrossover%2==1) numForCrossover++;
        for(int i=0;i<numForCrossover/2;i++){
            boolean notFound=true;
            int index = 0;
            while(notFound){
                index= (int) (Math.random()*(this.initialPopulationSize));
                if(!indicesSelected.contains(index)){
                    notFound=false;
                    indicesSelected.add(index);
                    list.add(index);
                }
            }
        }
        return list;
    }
    private void performCrossover(Map<Integer,Integer> map){
        for(Integer i:map.keySet()){
            double d=Math.random();
            if(d>this.crossoverProbability){
                continue;
            }
            int otherIndex=map.get(i);
            int a=fitness[i];
            int b=fitness[otherIndex];
            int maxIndex=i;
            int maxFit=a;
            if(b>a) {
                maxIndex=otherIndex;
                maxFit=fitness[maxIndex];
            }
            char[] c=this.individuals[maxIndex];
            this.cross(i,otherIndex);
            //System.out.println("hit line 451");
            double db=Math.random();
            if(db>0.5){
                int worstIndex=this.maxFit.maxIndex();
                int fitnessOfWorstIndex=this.fitness[worstIndex];
                this.individuals[worstIndex]=c;
                this.fitness[worstIndex]=maxFit;
                this.totalFitness=this.totalFitness-fitnessOfWorstIndex+maxFit;
                this.maxFit.changeKey(worstIndex,maxFit);
            }
        }
    }
    private void cross(int one, int two) { //finish figuring out how to mate them, maybe first 1/2 of alphabet
        char[] a=this.individuals[one];
        char[] b=this.individuals[two];
      // System.out.println("\n\n crossing index: "+one+" with index: "+two);
        int fitOldA=fitness[one];
        int fitOldB=fitness[two];
        char[] newA=new char[10];
        char[] newB=new char[10];
        List<Character> chars=new ArrayList<>(this.uniqueChars);
        int size=chars.size()/2;
        for(int i=0;i<size;i++){
            char c=chars.get(0);
            newA[this.find(a,c)]=c;
            newB[this.find(b,c)]=c;
            chars.remove(0);
        }
        while(chars.size()>0){
            char c=chars.remove(0);
            int indexForNewA=this.find(b,c);
            int indexForNewB=this.find(a,c);
            int current=indexForNewA;
            boolean placed=false;
            while(!placed){
                if(newA[current]=='\u0000'){
                    newA[current]=c;
                    placed=true;
                }
                else{
                    current++;
                    if(current>newA.length-1) current=0;
                }
            }
            current=indexForNewB;
            placed=false;
            while(!placed){
                if(newB[current]=='\u0000'){
                    newB[current]=c;
                    placed=true;
                }
                else{
                    current++;
                    if(current>newB.length-1) current=0;
                }
            }
        }
        this.individuals[one]=newA;
        this.individuals[two]=newB;
        int fitNewA=this.getFitness(one);
        int fitNewB=this.getFitness(two);
        this.fitness[one]=fitNewA;
        this.fitness[two]=fitNewB;
        this.totalFitness=this.totalFitness-fitOldA-fitOldB+fitNewA+fitNewB;
        this.maxFit.changeKey(two,fitNewB);
        if(fitNewA==0){
            this.foundSolution=true;
            this.solutionIndex=one;
        }
        else if(fitNewB==0){
            this.foundSolution=true;
            this.solutionIndex=two;
        }
    }
    private int find(char[] chars,char c){
        for(int i=0;i<chars.length;i++){
            if(chars[i]==c) return i;
        }
        return -1;
    }
    private Map<Integer,Integer> rouletteSelectionForCrossover(){
        Map<Integer,Integer>map=new HashMap<>();
        Set<Integer> indicesSelected=new HashSet<>();
        int currentIndex=0;
        int[] wheel=new int[1000];
        for(int k=0;k<1000;k++){
            wheel[k]=Integer.MAX_VALUE;
        }
        for(int i=0;i<this.fitness.length;i++){
            int fit=this.fitness[i];
            double percentAsDouble=(this.totalFitness/fit)*10.0; //to better deal with rlly low with percentages
            int percentAsInt=(int)percentAsDouble;
            for(int j=currentIndex;j<currentIndex+percentAsInt;j++){
                wheel[j]=i;
            }
            currentIndex+=percentAsInt;
        }
        int numForCrossover=this.initialPopulationSize/10;
        if(numForCrossover%2==1) numForCrossover++;
        for(int i=0;i<numForCrossover/2;i++) {
            boolean notFound = true;
            int index1 = 0, index2 = 0;
            int indexLeft=0,indexRight=0;
            while (notFound) {
                index1 = (int) (Math.random() * (1000));
                if (!indicesSelected.contains(wheel[index1])) {
                    notFound = false;
                    indicesSelected.add(wheel[index1]);
                    indexLeft=wheel[index1];
                }
            }
            notFound=true;
            while(notFound){
                index2= (int) (Math.random()*(1000));
                if(!indicesSelected.contains(wheel[index2])){
                    notFound=false;
                    indicesSelected.add(wheel[index2]);
                    indexRight=wheel[index2];
                }
            }
            map.put(indexLeft,indexRight);
        }
        return map;
    }
    private Map<Integer,Integer> tournamentSelectionForCrossover(){
        Map<Integer,Integer>map=new HashMap<>();
        Set<Integer> indicesSelected=new HashSet<>();
        int numForCrossover=this.initialPopulationSize/10;
        if(numForCrossover%2==1) numForCrossover++;
        for(int i=0;i<numForCrossover/2;i++){
            boolean notFound=true;
            int index1 = 0,index2 =0,index3=0,index4 =0;
            while(notFound){
                index1= (int) (Math.random()*(this.initialPopulationSize));
                if(!indicesSelected.contains(index1)){
                    notFound=false;
                    indicesSelected.add(index1);
                }
            }
            notFound=true;
            while(notFound){
                index2= (int) (Math.random()*(this.initialPopulationSize));
                if(!indicesSelected.contains(index2)){
                    notFound=false;
                    indicesSelected.add(index2);
                }
            }
            notFound=true;
            while(notFound){
                index3= (int) (Math.random()*(this.initialPopulationSize));
                if(!indicesSelected.contains(index3)){
                    notFound=false;
                    indicesSelected.add(index3);
                }
            }
            notFound=true;
            while(notFound){
                index4= (int) (Math.random()*(this.initialPopulationSize));
                if(!indicesSelected.contains(index4)){
                    notFound=false;
                    indicesSelected.add(index4);
                }
            }
            int indexLeft,indexRight;
            if(this.fitness[index1]>=this.fitness[index2]) indexLeft=index1;
            else indexLeft=index2;
            if(this.fitness[index3]>=this.fitness[index4]) indexRight=index3;
            else indexRight=index4;
            map.put(indexLeft,indexRight);
        }
        return map;
    }
    private int getFitness(int index){
        char[] chars=this.individuals[index];
        Map<Character,Integer> map=new HashMap<>();
        for(int i=0;i<10;i++){
            map.put(chars[i],i);
        }
        String au="";
        String ad="";
        String s="";
        for(int j=0;j<this.augend.length();j++){
            au+=map.get(this.augend.charAt(j));
        }
        for(int x=0;x<this.addend.length();x++){
            ad+=map.get(this.addend.charAt(x));
        }
        for(int y=0;y<this.sum.length();y++){
            s+=map.get(this.sum.charAt(y));
        }
        return Math.abs(Integer.parseInt(s)-Integer.parseInt(au)-Integer.parseInt(ad));
    }
    private void initializeFitnesses(){
        for(int i=0;i<this.initialPopulationSize;i++){
            int fitness=this.getFitness(i);
            this.fitness[i]=fitness;
            this.totalFitness+=fitness;
            this.maxFit.insert(i,fitness);
            if(fitness==0){
                this.foundSolution=true;
                this.solutionIndex=i;
                break;
            }

        }
    }
    private void initializeVars(GeneticAlgorithmConfig gac){
        this.crossoverProbability=gac.getCrossoverProbability();
        this.mutationProbability=gac.getMutationProbability();
        this.initialPopulationSize=gac.getInitialPopulationSize();
        this.maxGenerations=gac.getMaxGenerations();
        this.maxFit=new IndexMaxPQ<>(this.initialPopulationSize);
        this.fitness=new int[this.initialPopulationSize];
        for(int i=0;i<this.fitness.length;i++){
            fitness[i]=Integer.MAX_VALUE;
        }
        this.individuals= new char[this.initialPopulationSize][];
        for(int i=0;i<this.initialPopulationSize;i++){
            this.individuals[i]=new char[10];
        }
    }
    public void createRandomPopulation(){
        for(int i=0;i<this.initialPopulationSize;i++){
            Set<Character> chars=new HashSet<>(this.uniqueChars);
            char[] individual=new char[10];
            for(Character c:chars){
                boolean found=false;
                int index=0;
                while(!found){
                    index=(int)(Math.random()*10);
                    if(individual[index]=='\u0000'){
                        found=true;
                    }
                }
                //System.out.println("681-out of while");
                individual[index]=c;
            }
            this.individuals[i]=individual;
        }
        //this.print();
    }
}
