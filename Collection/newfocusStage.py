import serial
import time

class smc100:
    def __init__(self,comPort):
        self.ser = serial.Serial(comPort, 57600, bytesize = 8, stopbits=1,parity = 'N', xonxoff=True,timeout=0.050)


    def get_position(self):
        command = '1TP?\r\n'.encode()
        self.ser.flush()
        self.ser.write(command)
        return float(str(self.ser.readline()).split('TP')[1].split('\\r')[0])
    

    def is_moving(self):
        #returns true if moving
        command = '1TS?\r\n'.encode()
        self.ser.write(command)
        try:
            output = str(self.ser.readline()).split('TS')[1].split('\\r')[0]
        except:
            return True
        return output == "000028"
        
        

    def move_absolute(self,position):
        #units of mm
        command= str('1PA' + str(position) + '\r\n').encode()
        self.ser.write(command)
        
    def wait_till_done(self):
        while self.is_moving():
            time.sleep(.05)
    
    def close(self):
        self.ser.close()
        

if __name__ == "__main__":
    t = smc100("COM1")
    t.move_absolute(12.4) #13.5 best
    while t.is_moving():
        print(str(t.get_position()))
        time.sleep(2)
        
    t.close()
    
