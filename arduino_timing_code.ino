/*
@author: Stanisław Niziński
################ THIS IS TIMING SCRIPT FOR ATMEL MEGA32U4        ##############
################ MICROCONTROLLERS (PREFERABLY ON ARDUNIO BOARD)  ##############
################ It is designed to generate timing sequence for  ##############
################ the multi-pulse transient absorption experiment ##############
################ with precision determined by MCU clock (usually ##############
################ 16MHz). Sequence must be initialized by RS232   ##############
################ message. Performance of this script is          ##############
################ tested, but it is also work in progress and     ##############
################ tailored for this specific experiment.          ##############
*/

//Firstly declare some working variables. In general, the idea is that hardware counter should run with low prescaller rate, so high time resolution
//is achieved. There is also second counter - software oriented one. Every time hardware-counter reaches the end, short interruption function is
//executed and software counter is updated. It can also update counter registers on the fly. Because of that, pulse edge must be generated after 
//some number of the clock cycles, because counter run continuously during interruption script, and script must be completed before any new counter-
//related event occurs. This procedure and external python script that uses it is designed such that counters are operated properly. Therefore, in order
//to use this script in a new way, one should ensure that counters are still happy and also carefully check timing using oscilloscope!

//for counters
volatile long overflow_count; //numbers of subsequent overflows incremented each time
volatile long overflow_max; //how many counter overflows has to be completed

//to define when pulses should be produced. there are two outputs, one for camera trigger one for laser trigger.
volatile long * laser_pulses_no = NULL; //array with numbers of overflows when pulse should be produced. first counter run is zero, but should not be used!
volatile long * camera_pulses_no = NULL; //array with numbers of overflows when pulse should be produced. first counter run is zero, but should not be used!
unsigned laser_pulse_count = 0; //counter, so that every time compare between given overflow number and array above is true, one can move to next number
unsigned camera_pulse_count = 0; //counter, so that every time compare between given overflow number and array above is true, one can move to next number
unsigned laser_pulse_size = 0; //number of numbers of pulses in array
unsigned camera_pulse_size = 0; //number of numbers of pulses in array

volatile long led_on_no = -1; //here possible to set led on at specific overflow interruption to avoid talking to arduino during counting
volatile long led_off_no = -1; //the same but to switch it off, in both these -1 means it is not in use

volatile bool pulsetrain = false; //if true, then at every overflow there will be a pulse!

//for uart
const unsigned max_chars = 128;
volatile unsigned last_char = 0;
volatile char buffer[max_chars];

void setupTwoPulsePWM (float period, float laser_pulse, float camera_pulse, long total_overflows) {
  //this function sets up counters and initializes all crucial registers
  //period - all cycles time - pulses spacing - in seconds
  //I tested periods down to 512 us - below stuff may be not fast enough...
  //laser_pulse, camera_pulse - HIGH pulse length
  //camera LOW must be 2-3 cycle or more! cannot be zero! otherwise first overflows wont work! probably 4 cycles (64us) is safe minimum!
  //total_overflows is total number of overflows before counter is stopped

  if(laser_pulses_no == NULL || laser_pulse_size == 0) {
    Serial.println("Error! Arrays with laser pulse numbers not allocated or empty!");
    return;
  }
  if(camera_pulses_no == NULL || camera_pulse_size == 0) {
    Serial.println("Error! Arrays with camera pulse numbers not allocated or empty!");
    return;
  }

  Serial.print("TCNT1: ");
  Serial.print(TCNT1);
  Serial.print("\n");

  laser_pulse_count == 0;
  camera_pulse_count == 0;

  overflow_count = 0;
  overflow_max = total_overflows;

  TCCR1A = 0; // reset timer 1 control register          
  TCCR1B = 0;

  GTCCR |= (1 <<PSRASY); // reset prescaler

  TCCR1B |= (1 << WGM12); // set TOP = ICR1 (WGM13:0 = 14 = 1110)
  TCCR1B |= (1 << WGM13);
  TCCR1A |= (1 << WGM11);

  //DDRB |= (1 << 1); // set OC1A as output! old atmega
  DDRB |= (1 << 5); // set OC1A (laser) as output! 32u4
  DDRB |= (1 << 6); // set OC1B (camera) as output too! 32u4
  
  //float prescaler = 1024.0; //modify TCCR1B to change this, currently it gives about 4s window and 64us resolution
  //float prescaler = 256.0; //modify TCCR1B to change this, currently it gives about 1048ms window and 16us resolution
  float prescaler = 64.0; //modify TCCR1B to change this, currently it gives about 262ms window and 4us resolution

  ICR1 = (unsigned int)( (16000000.0f*period)/prescaler-1.0f ); //set TOP value
  
  Serial.print("ICR1: ");
  Serial.print(ICR1);
  Serial.print("\n");

  //print exact length of one period - start to overflow
  float exact_period = prescaler*(1+float(ICR1))/16000000.0f;
  Serial.print("exact period: ");
  Serial.print(exact_period,9);
  Serial.print("\n");

  //print exact length of one cycle
  float exact_cycle = prescaler/16000000.0f;
  Serial.print("exact cycle: ");
  Serial.print(exact_cycle,9);
  Serial.print("\n");

  OCR1A = ICR1-(unsigned int)( laser_pulse/(prescaler/16000000.0f) ); //set compare math value for PWM A (where HIGH starts) - laser

  //low will last 1+OCR1A cycles
  //high will last then ICR1+1 - (OCR1A+1) cycles

  Serial.print("OCR1A: ");
  Serial.print(OCR1A);
  Serial.print("\n");

  //print exact LOW time of laser channel in last 2 overflows
  float exact_laser_low = prescaler*(1+float(OCR1A))/16000000.0f;
  Serial.print("exact laser LOW: ");
  Serial.print(exact_laser_low,9);
  Serial.print("\n");

  OCR1B = ICR1-(unsigned int)( camera_pulse/(prescaler/16000000.0f) ); //set compare math value for PWM A (where HIGH starts) - camera

  //print exact LOW time of camera channel
  float exact_camera_low = prescaler*(1+float(OCR1B))/16000000.0f;
  Serial.print("exact camera LOW: ");
  Serial.print(exact_camera_low,9);
  Serial.print("\n");

  TCNT1 = 0; // clock value to zero!
  TIMSK1 |= (1 << ICIE1);   // interrupt on Timer 1 input capture (ICT match == TOP)
  //TIMSK1 |= (1 << OCIE1A);  //diagnostics - on compare match


  //TCCR1B |= (1 << CS10) | (1 << CS12); //prescaller 1024 and start timer! //every 64us
  //TCCR1B |= (1 << CS12); //prescaller 256 and start timer! //every 16us
  TCCR1B |= (1 << CS10) | (1 << CS11); //prescaller 64 and start timer! //every 4us
  
}

//ISR (TIMER1_COMPA_vect) {  //on compare mathch - eventually for diagnostics
//  Serial.println(TCNT1);
//  Serial.println(OCR1A);
//}

ISR (TIMER1_CAPT_vect) { //timer capture interruption (works as overflow, but must change then ICIE1 to TOIE1)
  overflow_count++;     // count number of Counter1 overflows  

  //firstly clear all pulse flags by default, and enable only appropriate ones for given overflow
  //note that counters are running continuously - this function must be FAST and counter compare match cannot occur too fast!
  if(pulsetrain == false){ 
    TCCR1A &= ~(1 << COM1A1);
    TCCR1A &= ~(1 << COM1A0);
    TCCR1A &= ~(1 << COM1B1);
    TCCR1A &= ~(1 << COM1B0);
  }

  if(overflow_count >= overflow_max){ //stop timer after given number of pulses
    TCCR1A = 0;
    TCCR1B = 0;
    Serial.println("seq done!");
    Serial.println(overflow_count);
    laser_pulse_count = 0;
    camera_pulse_count = 0;
  }

  else if(overflow_count == laser_pulses_no[laser_pulse_count]){ //means that next one should produce laser pulse
    TCCR1A |= (1 << COM1A1); // OC1A (laser) will go HIGH at compare match! and LOW at bottom! (PWM)
    TCCR1A |= (1 << COM1A0);

    if(laser_pulse_count < laser_pulse_size-1) laser_pulse_count++;

  }

  else if(overflow_count == camera_pulses_no[camera_pulse_count]){ //means that next one should produce camera pulse
    TCCR1A |= (1 << COM1B1); // OC1B (camera) will go HIGH at compare match! and LOW at bottom! (PWM)
    TCCR1A |= (1 << COM1B0); 

    if(camera_pulse_count < camera_pulse_size-1) camera_pulse_count++;

  }
  
  //at this point of function execution, counter will probably make a few counts. one must be sure that compare match is 
  //scheduled later, not before (so that it is not skipped)!

  if(overflow_count == led_on_no){ //can be used to switch led on
    PORTD |= (1 << 4);
  }
  if(overflow_count == led_off_no){ //can be used to switch led off
    PORTD &= ~(1 << 4);
  }
      
}



void setup() { //just setup UART and outputs
  DDRD |= (1 << 4); // set PORTD4 as output for led triggers
  PORTD &= ~(1 << 4);
  
  Serial.begin(9600);
  Serial.println("Initialized!");
}

void loop() {
  //capture uart - look for commands and do what has to be done...
  char rc;

  if (Serial.available() > 0) {
    rc = Serial.read();

    if (rc != '\n' && last_char <= max_chars-1) { //save char
      buffer[last_char] = rc;
      last_char++;
      
    }
    else if (rc == '\n' && last_char <= max_chars-1) { //end of line! process it!
      buffer[last_char] = '\0'; //terminate string
  
      char buffer_tmp[max_chars];
      strcpy(buffer_tmp, buffer); //String does not want to eat volatile char, I had to use workaround.
      String input = String(buffer_tmp);
      int nexpected = 25;
      String substrings[nexpected];
      
      int tab_index = 0;
      int index = 0;
      
      //process max of nexpected space-separated words in string
      while (index < input.length()) {
        int spaceIndex = input.indexOf(' ', index); //find space
        String substring;

        if (spaceIndex == -1) { // If no more spaces, use the rest of the string
          substring = input.substring(index);
        }
        else { // Extract the substring between the current index and the space
          substring = input.substring(index, spaceIndex);
        }
        
        
        if(tab_index < nexpected) { //save substring if there is place to do it
          substrings[tab_index] = substring;
          tab_index++;
        }

        // Move to the next substring or break
        if (spaceIndex == -1) break;
        else index = spaceIndex + 1;
      }

      
      if (tab_index > 1) { //setup sequence for pulses
        int command = substrings[0].toInt();

        //CHECK FOR COMMANDS BELOW

        //RUN SEQUENCE COMMAND ----------------------------------------------------
        if(command == 3 && tab_index == 5) {
          float period = substrings[1].toFloat();
          float laser_pulse = substrings[2].toFloat();
          float camera_pulse = substrings[3].toFloat();
          long total_overflows = substrings[4].toInt();
		  //in this command, one specifies i) period - how long one timer run will take,
		  //ii) laser_pulse - length of laser trigger, iii) length of the camera trigger, 
		  //iv) total_overflows - how many runs of the timer will occur in one sequence,
		  //so whole sequence will take periods x total_overflows
		  //Note, that when timer begins counting, output is LOW, and after compare match,
		  //it is HIGH until the end. This way is better (so interruption func has time to update regs)

          Serial.println("Starting sequence!");
          Serial.print(period, 9);
          Serial.print("\n");
          Serial.print(laser_pulse, 9);
          Serial.print("\n");
          Serial.print(camera_pulse, 9);
          Serial.print("\n");
          Serial.print(total_overflows);
          Serial.print("\n");
          setupTwoPulsePWM(period, laser_pulse, camera_pulse, total_overflows);
        }

        //SET LASER PULSE NUMBERS COMMAND ----------------------------------------------------
        else if(command == 1 && tab_index > 2) {
			//this is just another command, to be sent before run sequence command (3). it defines
			//specific counter runs for which trigger has to be generated. This way we can keep
			//high timing precision, and long durations between pulses. So just for some counter
			//overloads, no trigger pulse will be generated. Note also that one cannot schedule too
			//many pulses, because memory is limited (but tens of them should be no problem).
          laser_pulse_size = substrings[1].toInt();

          if(laser_pulses_no != NULL) {
            free(laser_pulses_no);
          }

          laser_pulses_no = (long*) malloc(laser_pulse_size * sizeof(long));

          if(laser_pulses_no == NULL) {
            laser_pulse_size = 0;
            Serial.println("Error! Cannot allocate array!");
          }
          else{
            for(unsigned i = 0; i < laser_pulse_size; i++){
              if(i+2 >= tab_index) {
                Serial.println("Error! Not enough space for numbers! Sequence will be invalid!");
                break;
              }
              laser_pulses_no[i] = substrings[i+2].toInt();
            }
          }
          Serial.println("Laser pulse seq. updated!");

        }

        //SET CAMERA PULSE NUMBERS COMMAND ----------------------------------------------------
        else if(command == 2 && tab_index > 2) {
			//this is the same as command 1 but sets pulses for camera.
          camera_pulse_size = substrings[1].toInt();

          if(camera_pulses_no != NULL) {
            free(camera_pulses_no);
          }

          camera_pulses_no = (long*) malloc(camera_pulse_size * sizeof(long));

          if(camera_pulses_no == NULL) {
            camera_pulse_size = 0;
            Serial.println("Error! Cannot allocate array!");
          }
          else{
            for(unsigned i = 0; i < camera_pulse_size; i++){
              if(i+2 >= tab_index) {
                Serial.println("Error! Not enough space for numbers! Sequence will be invalid!");
                break;
              }
              camera_pulses_no[i] = substrings[i+2].toInt();
            }
          }
          Serial.println("Camera pulse seq. updated!");
          
        }

        //LED COMMAND ----------------------------------------------------
        else if(command == 0 && tab_index == 2) { //do stuff to LED output
		  //this is just to manually start/stop LED at command execution

          int led_on = substrings[1].toInt();

          if(led_on == 1) {
            PORTD |= (1 << 4);
            Serial.println("LED is now on!");
          }
          else if(led_on == 0) {
            PORTD &= ~(1 << 4);
            Serial.println("LED is now off!");
          }

        }

        //AUTO LED ON OFF COMMAND ----------------------------------------------------
        else if(command == 4 && tab_index == 3) { //if led should be changed during counting, use this instead
		  //this is just additional feature, to start led at some specific counter run
		  //this way does not have precise accuracy, but it is not needed for LED, it
		  //just needs to be started before camera starts running...

          led_on_no = substrings[1].toInt(); //first number is overflow number when to led on
          led_off_no = substrings[2].toInt(); //second is when to led off. in both cases, -1 disables this feature

          Serial.println("Updated LED at overflow counts!");

        }

        //PULSE TRAIN OVERRIDE ----------------------------------------------------
        //just for pulse power measurement/diagnostics - runs pulses at all oveflows, ignores commands 1 and 2
        else if(command == 5 && tab_index == 2) {
          int override_on = substrings[1].toInt();

          if(override_on == 1) pulsetrain = true;
          else if(override_on == 0) pulsetrain = false;

          Serial.println("Laser override updated!");

        }

        else{
		  //not recognized command, something went wrong!
          Serial.println("Serial error! Wrong command!");
        }

        //COMMANDS BLOCK ENDS HERE

      }

      else {
		//not recognized command, something went wrong!
        Serial.println("Serial error! Wrong command!");
      }

      last_char = 0;
    }
    else {
	  //error! too long incoming line!
      last_char = 0;
      Serial.println("Serial error!");
    }

  }

}


