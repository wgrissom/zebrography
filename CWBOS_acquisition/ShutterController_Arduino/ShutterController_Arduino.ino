int doi;
float time_delay;
void setup() {
  // put your setup code here, to run once:
  pinMode(A1,OUTPUT);
  pinMode(A2,OUTPUT);
  pinMode(A3,OUTPUT);
  pinMode(10,OUTPUT);
  Serial.begin(9600);
}

void loop() {
  // put your main code here, to run repeatedly:
  digitalWrite(A1,LOW);
  digitalWrite(A2,LOW);
  digitalWrite(A3,LOW);
  digitalWrite(10,LOW);
  delay(100);
  if (Serial.available()>0)
  {
    doi = Serial.read();
    doi = doi - 48;
    delay(200);
    time_delay = Serial.parseFloat(); // wait for time_delay ms 
//    Serial.println(doi);
    delay(200);
    if (doi == 0)
    {
      digitalWrite(A1,LOW);
      digitalWrite(A2,LOW);
      digitalWrite(A3,LOW);
      digitalWrite(10,LOW); // do not trigger the waveform generator
      
      digitalWrite(A1,HIGH);
      digitalWrite(A2,HIGH);
      digitalWrite(A3,HIGH); //turn the camera by set A1,A2 and A3 as "high"
      delay(400); // wait for camera to take a pictures. 
      digitalWrite(A1,LOW);
      digitalWrite(A2,LOW);
      digitalWrite(A3,LOW);

    }
    else if (doi == 1)
    {
      digitalWrite(A1,LOW);
      digitalWrite(A2,LOW);
      digitalWrite(A3,LOW);
      digitalWrite(10,LOW);
      digitalWrite(A1,HIGH); // 
      digitalWrite(A2,HIGH); //
      digitalWrite(A3,HIGH); // turn on the camera
      delay(time_delay); // wait for several miliseconds and then trigger the waveform generator. 
      digitalWrite(10,HIGH);
      delay(400);
      digitalWrite(A1,LOW);
      digitalWrite(A2,LOW);
      digitalWrite(A3,LOW);
      digitalWrite(10,LOW);
    }
  }
}
